#include <array>
#include <math.h>
#include "ascent/Ascent.h"

const double PWM_FREQUENCY = 144e6; // 36kHz update rate
const unsigned MAX_PERIOD = 2048;

struct d_eps {
    static constexpr double eps = 1e-9;
};

using Sampler = asc::SamplerT<double, d_eps>;

using namespace std;
struct CAPWM {
    double tq;
    double tp;
    double ts;
    bool up;
    unsigned period;
    array<unsigned, 4> chan;
    double t0;
    array<double, 3> Vabc;
    CAPWM(double fclk, unsigned p) : up(false), period(p), chan{0},  t0(0.0) {
        tq = 1.0/fclk;
        tp = period * tq;
        ts = 2 * tp;
        chan[3] = period - 1;
        Vabc[0] = Vabc[1] = Vabc[2] = 0.0;
    }
    void setSamplePoint(unsigned x) {chan[3] = x;}
    bool operator()(Sampler &sampler, double &t) {
        if (sampler(tp)) {
            up = !up;
            t0 = t;
        }
        if (up) {
            for (unsigned i = 0; i < 3; ++i)
                if (sampler.event(t0 + (period - chan[i]) * tq))
                    Vabc[i] = 1.0;
            return false;
        } else {
            for (unsigned i = 0; i < 3; ++i)
                if (sampler.event(t0 + chan[i] * tq))
                    Vabc[i] = 0;
            return (sampler.event(t0 + chan[3] * tq));
        }
    }
};

template <typename T>
static inline tuple<T, T, T> iclarke(const T a, const T b)
{
    return make_tuple(a, (sqrt(3.0)*b - a) * 0.5, -(sqrt(3.0)*b + a) * 0.5);
}

template <typename T>
static inline void ipark(const double sin_phi, const double cos_phi, const array<T, 2> &dq, array<T, 2> &ab)
{
    ab[0] = dq[0] * cos_phi - dq[1] * sin_phi;
    ab[1] = dq[1] * cos_phi + dq[0] * sin_phi;
}

template <typename T>
struct PIRegulator {
    T Ka, Kb;
    T limit;
    T sum;
    PIRegulator(): Ka{0}, Kb{0}, sum(0) {}
    T operator()(const T err, const T ff = 0.0) {
        T out = Ka * err;
        T x = Kb * out;
        out += ff;
        T mx = fmax(limit - out, 0);
        T mn = fmin(-limit - out, 0);
        sum += x;
        if (sum > mx)
            sum = mx;
        if (sum < mn)
            sum = mn;
        out += sum;
        if (out > limit)
            out = limit;
        if (out < -limit)
            out = -limit;
        return out;
    }
};

template <typename T>
struct PIRegS {
    T Ka, Kb;
    T sum;
    PIRegS(): Ka{0}, Kb{0}, sum(0) {}
    T operator()(const T err) {
        T out = Ka * err;
        sum += Kb * out;
        out += sum;
        return out;
    }
};

// RP34M-221V48 540W NEMA34 motor
static const double VBUS = 48;
static const double Rs = 0.15;
static const double Ls = 3.4e-4;
static const double Ke = 0.074;
static const double J = 2.74e-4;
static const double Ja = 0.0;
static const double Fv = 0.0;
static const double Kf = 5e-6;
static const double poles = 4;

template <typename T>
struct PIReg {
    T Kp, Ki;
    T sum;
    PIReg(): Kp{0}, Ki{0}, sum(0) {}
    T operator()(const T err) {
        sum += Ki * err;
        return sum + Kp * err;
    }
};

template <typename T>
static inline void clarke3(const T a, const T b, const T c,  T &al, T &bt)
{
    al = (a * 2.0 - (b + c)) / 3;
    bt = (b - c)/sqrt(3);
}

template <typename T>
static inline void clarke3(const array<T, 3> &abc, array<T, 2> &ab)
{
    clarke3(abc[0], abc[1], abc[2], ab[0], ab[1]);
}

template <typename T>
static inline void park(const double sin_phi, const double cos_phi, const array<T, 2> &ab, array<T, 2> &dq)
{
   dq[0] = ab[0] * cos_phi + ab[1] * sin_phi;
   dq[1] = ab[1] * cos_phi - ab[0] * sin_phi;
}

const double R1 = 82e3;
const double R2 = 5e3;
const double C = 22e-9;

const double VSCALE_DIV = (R1 + R2) / R2;

struct Motor {
    asc::Param Id, Iq, oe, theta;
    double Vbus;
    array<double, 3> Vabc;
    double Ia, Ib, Ic;
    Motor(asc::state_t& state) : Id(state), Iq(state),
        oe(state), theta(state) {
        Iq = Id = 0.0;
        oe = 0.0;
        theta = 0.0;
    }
    void operator()(const asc::state_t&, asc::state_t& D, const double) {
        array<double, 2> Vab;
        array<double, 2> Vdq;
        clarke3(Vabc, Vab);
        park(sin(theta), cos(theta), Vab, Vdq);
        Id(D) = (Vdq[0] + oe * Ls * Iq - Rs * Id) / Ls;
        Iq(D) = (Vdq[1] - oe * (Ls * Id + Ke) - Rs * Iq) / Ls;

        double Te = 1.5 * Ke * poles * Iq;
        double Tf = oe;
        Tf = (Tf > Kf) ?  Kf : Tf;
        Tf = (Tf < -Kf)? -Kf : Tf;
        oe(D) = poles * (Te - Tf - Fv * oe / poles) / (J + Ja);
        theta(D) = oe;

        array<double, 2> Idq;
        array<double, 2> Iab;
        Idq[0] = Id;
        Idq[1] = Iq;
        ipark(sin(theta), cos(theta), Idq, Iab);
        tie(Ia, Ib, Ic) = iclarke(Iab[0], Iab[1]);
    }
};


struct Hall120 {
    enum HAL_VAL {
        HAL_0_60 = 0b100,
        HAL_60_120 = 0b110,
        HAL_120_180 = 0b010,
        HAL_180_240 = 0b011,
        HAL_240_300 = 0b001,
        HAL_300_360 = 0b101
    };
    HAL_VAL val;
    Hall120(): val(HAL_0_60) {}
    bool operator()(double angle) {
        HAL_VAL v = HAL_0_60;
        if (angle < M_PI/3)
            v = HAL_0_60;
        else if (angle < 2*M_PI/3)
            v = HAL_60_120;
        else if (angle < 3*M_PI/3)
            v = HAL_120_180;
        else if (angle < 4*M_PI/3)
            v = HAL_180_240;
        else if (angle < 5*M_PI/3)
            v = HAL_240_300;
        else
            v = HAL_300_360;
        bool change = v != val;
        val = v;
        return change;
    }
};

struct Simulator;
class System {
    CAPWM pwm;
    Motor motor;
    Hall120 hes;
    friend Simulator;
public:
    asc::Param Vca, Vcb, Vcc;
    System(asc::state_t &s) :  pwm(PWM_FREQUENCY, MAX_PERIOD), motor(s), Vca(s), Vcb(s), Vcc(s) {
        motor.Vbus = VBUS;
        Vca = Vcb = Vcc = 0.0;
    }
    double getTheta() {return motor.theta;}
    double getOmega() {return motor.oe;}
    void getCurrents(array<double, 3> &iabc) {iabc = {motor.Ia, motor.Ib, motor.Ic};}
    double getTs() {return pwm.ts;}
    void setPWM(array<unsigned, 3> &x) {
        for (unsigned i = 0; i < 3; ++i)
            pwm.chan[i] = x[i];
    }
    void operator()(const asc::state_t& x, asc::state_t& D, const double t) {
        for (unsigned i = 0; i < 3; ++i)
            motor.Vabc[i] = motor.Vbus * pwm.Vabc[i];
        motor(x, D, t);
        Vca(D) = ((motor.Vabc[0] - Vca) / R1 - Vca / R2) / C;
        Vcb(D) = ((motor.Vabc[1] - Vcb) / R1 - Vcb / R2) / C;
        Vcc(D) = ((motor.Vabc[2] - Vcc) / R1 - Vcc / R2) / C;
    }
};

struct fo_filter
{
    float a, b0, b1, y1, x1;
    fo_filter(float beta) {
        a = (beta - 2.0) / (beta + 2.0);
        b1 = b0 = beta / (beta + 2.0);
        y1 = x1 = 0;
    }
    float operator() (const float inputValue) {
        float y0 = (b0 * inputValue) + (b1 * x1) - (a * y1);
        x1 = inputValue;
        y1 = y0;
        return y0;
    }
};

template <typename T>
struct AngleEstimator {
    const T Ts;
    T Vr_a, Vr_b;
    T Va, Vb;
    T Vemf_a, Vemf_b;
    T angle;
    T omega;
    T angle_err;
    T bemf_square;
    fo_filter omega_f;
    PIRegS<T> omega_pi;
    PIRegS<T> pi_bemf_a, pi_bemf_b;
    AngleEstimator(const T _Ts) : Ts(_Ts), omega_f(0.002) {
        Vr_a = Vr_b = Vemf_a = Vemf_b = 0;
        omega = angle = 0.0;
        pi_bemf_a.Ka = pi_bemf_b.Ka = Ls * 4000;
        pi_bemf_a.Kb = pi_bemf_b.Kb = Ts*Rs/Ls;
        omega_pi.Ka = 100;
        omega_pi.Kb = Ts*40;
    }
    bool operator()(const array<double, 2> Vab, const array<double, 2> Iab){
        const double k = (Rs / Ls) * Ts;
        Vr_a += k * (Vab[0] - Vemf_a - Vr_a);
        Vr_b += k * (Vab[1] - Vemf_b - Vr_b);
        Vemf_a = pi_bemf_a(-Iab[0] + Vr_a/Rs);
        Vemf_b = pi_bemf_b(-Iab[1] + Vr_b/Rs);
        angle_err = -2.0*Vemf_a*Vemf_b*cos(2*angle)+
                (Vemf_a*Vemf_a-Vemf_b*Vemf_b)*sin(2*angle);
        bemf_square = Vemf_a*Vemf_a + Vemf_b*Vemf_b;
        omega = omega_pi(angle_err);
        angle += omega * Ts;
        omega = omega_f(omega);
        while (angle > 2 * M_PI) angle -= 2 * M_PI;
        while (angle < 0) angle += 2 * M_PI;
        return (bemf_square > 0.1)?true:false;
    }
};

struct HallEstimator {
    const double Ts;
    double To;
    double theta_hall, theta_in;
    Hall120::HAL_VAL val;
    double angle;
    double omega;
    double err;
    PIRegulator<double> omega_pi;
    fo_filter omega_f;
    HallEstimator(const double ts) : Ts(ts), To(0), val(Hall120::HAL_0_60), omega_f(Ts*100) {
        theta_in = theta_hall = omega = angle = 0;
        omega_pi.Ka = 160;
        omega_pi.Kb = Ts * 20;
        omega_pi.limit = 2*M_PI*1000;
    }
    void operator()(const double t, Hall120::HAL_VAL v) {
        To = t;
        switch(v) {
        case Hall120::HAL_0_60:
            theta_hall = 0;
            break;
        case Hall120::HAL_60_120:
            theta_hall = M_PI/3;
            break;
        case Hall120::HAL_120_180:
            theta_hall = 2*M_PI/3;
            break;
        case Hall120::HAL_180_240:
            theta_hall = 3*M_PI/3;
            break;
        case Hall120::HAL_240_300:
            theta_hall = 4*M_PI/3;
            break;
        case Hall120::HAL_300_360:
            theta_hall = 5*M_PI/3;
            break;
        }
    }
    void operator()(const double t) {
        double dt = t - To;
        double dtheta = omega * dt;
        //if (dtheta > 2*M_PI/3) {
            dtheta = 0;
        //}
        theta_in = theta_hall + dtheta;
        /* sin(a-b) = cos(b)*sin(a) - cos(a)*sin(b) */
        err = -cos(theta_in)*sin(angle) + cos(angle)*sin(theta_in);
        omega = omega_pi(err);
        angle += omega * Ts;
        omega = omega_f(omega);
        while (angle > 2 * M_PI) {
            angle -= 2 * M_PI;
        }
        while (angle < 0) {
            angle += 2 * M_PI;
        }
    }
};

struct Control {
    array<double, 2> idq;
    double om_req;
    array<double, 2> vdq;
    array<double, 2> vab;
    array<double, 3> vabc;
    array<unsigned, 3> chan;
    PIRegulator<double> pi_id;
    PIRegulator<double> pi_iq;
    AngleEstimator<double> a_est;
    HallEstimator he;
    PIRegulator<double> speed_pi;
    enum FSM {
        IDLE = 0,
        RUN,
        REVERSE,
    } fsm;
    Control(double ts) : idq{0}, vdq{0}, vab{0}, vabc{0}, a_est (ts), he(ts) {
        pi_iq.Ka = pi_id.Ka = Ls * 2 * M_PI * 4000;
        pi_iq.Kb = pi_id.Kb = ts * Rs/Ls;
        pi_id.limit = pi_iq.limit = VBUS*1.15/2;
        speed_pi.Ka = 0.002;
        speed_pi.Kb = ts*.5;
        speed_pi.limit = 15.0;
        om_req = 0;
    }
    void svgen(const array<double, 2> &ab, array<double, 3> &abc, const double Vbus)
    {
        const double scale = 2.0 / Vbus;
        tie(abc[0], abc[1], abc[2]) = iclarke(ab[0], ab[1]);
        const double mn = fmin(abc[0], fmin(abc[1], abc[2]));
        const double mx = fmax(abc[0], fmax(abc[1], abc[2]));
        const double off = (mn + mx) * .5;
        for (unsigned i = 0; i < 3; ++i) {
            abc[i] -= off;
            abc[i] *= scale;
        }
    }
    void setDuty(const array<double, 3> d) {
        for (unsigned i = 0; i < 3; ++i) {
            double x = MAX_PERIOD * (d[i] + 1.0) * 0.5;
            if (x < 0)
                x = 0;
            if (x > MAX_PERIOD - 1)
                x = MAX_PERIOD - 1;
            chan[i] = int(x);
        }
    }
    void hall_irq(double t, Hall120::HAL_VAL v) {he(t, v);}
    void run(double t, System &sys) {
        array<double, 3> Iabc;
        sys.getCurrents(Iabc);
        array<double, 2> Iab, Vab, Idq;
        clarke3(Iabc, Iab);
        clarke3(sys.Vca*VSCALE_DIV, sys.Vcb*VSCALE_DIV, sys.Vcc*VSCALE_DIV, Vab[0], Vab[1]);
        //double a = sys.getTheta(), o = sys.getOmega();
        double a = he.angle, o = 0.0;
        if (t > 0.1)
            om_req = 200;
        if (t > 0.5)
            om_req = -200;
        if (t > 1.5)
            om_req = 0;
        bool est_valid = a_est(Vab, Iab);
        he(t);
        if (est_valid) {
            a = a_est.angle;
            o = a_est.omega;
        }
        park(sin(a), cos(a), Iab, Idq);
        idq[1] = speed_pi(om_req - o);
        // current loop
        vdq[0] = pi_id(idq[0] - Idq[0], -o * Idq[1] * Ls);
        vdq[1] = pi_iq(idq[1] - Idq[1], o * (Ls * Idq[0] + Ke));
        ipark(sin(a), cos(a), vdq, vab);
        tie(vabc[0], vabc[1], vabc[2]) = iclarke(vab[0], vab[1]);
        svgen(vab, vabc, VBUS);
        setDuty(vabc);
        sys.setPWM(chan);
    }
    void pwm_irq(double t, System &sys) {
        static unsigned fsm = 0;
        array<double, 3> Iabc;
        sys.getCurrents(Iabc);
        array<double, 2> Iab, Vab, Idq;
        clarke3(Iabc, Iab);
        clarke3(sys.Vca*VSCALE_DIV, sys.Vcb*VSCALE_DIV, sys.Vcc*VSCALE_DIV, Vab[0], Vab[1]);
        double a = he.angle, o = 0.0;
        vdq[0] = vdq[1] = 0.0;
        static double iprev = 0.0;
        static double ls = 0.0;
        if (t > 0.1) {
            if (fsm == 0) {
                vdq[0] = 0.2 * VBUS;
                fsm = 1;
            } else if (fsm == 1) {
                vdq[0] = 0;
                fsm = 2;
                iprev = Iab[0];
                ls = 0.2 * VBUS * he.Ts / Iab[0];
                printf("Ls is %f\n", ls);
            } else if (fsm == 2) {
                double di = iprev - Iab[0];
               static int i = 0;
                if (i++ > 8)
                    fsm = 3;
                double rs = ls * di/(i*he.Ts * iprev);
                printf("Rs is %f\n", rs);
            }
        }
        // current loop
        //vdq[0] = pi_id(idq[0] - Idq[0],-o * Idq[1] * Ls);
        //vdq[1] = pi_iq(idq[1] - Idq[1], o * (Ls * Idq[0] + Ke));
        ipark(sin(a), cos(a), vdq, vab);
        tie(vabc[0], vabc[1], vabc[2]) = iclarke(vab[0], vab[1]);
        svgen(vab, vabc, VBUS);
        setDuty(vabc);
        sys.setPWM(chan);
    }
};

struct Simulator {
    asc::state_t s;
    asc::RK4 integrator;
    asc::Recorder rec;
    System *sys;
    Simulator() {
        s.reserve(16);
        sys = new System(s);
    }
    void run(Control &ctrl, double t_end);
};

void Simulator::run(Control &ctrl, double t_end) {

    double t = 0.0, dt = 1e-4;
    ctrl.pwm_irq(t, *sys);
    //ctrl.om_req = 100;
    while (t < t_end) {
        Sampler sampler(t, dt);
        if (sys->pwm(sampler, t)) {
            ctrl.pwm_irq(t, *sys);
            rec({t, sys->getOmega(), sys->getTheta(), ctrl.a_est.omega, ctrl.a_est.angle,
                 ctrl.he.err, ctrl.he.theta_hall, ctrl.he.theta_in, ctrl.he.omega, ctrl.he.angle});
        }
        integrator(*sys, s, t, dt);
        while (sys->motor.theta > 2 * M_PI) {
            sys->motor.theta -= 2 * M_PI;
        }
        while (sys->motor.theta < 0) {
            sys->motor.theta += 2 * M_PI;
        }
        if (sys->hes(sys->motor.theta))
            ctrl.hall_irq(t, sys->hes.val);
    }
    rec.csv("pwm", {"t", "omega", "angle", "omega_est", "angle_est", "err_he",
                    "angle_hall", "angle_e", "omega_he", "hall_est"});
}

int main()
{
    Simulator sim;
    Control c(sim.sys->getTs());
    sim.run(c, 2);
    return 0;
}
