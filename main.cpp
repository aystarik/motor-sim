#include <array>
#include <math.h>
#include "ascent/Ascent.h"

const double PWM_FREQUENCY = 144e6; // 36kHz update rate
const unsigned MAX_PERIOD = 2048;

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
    void setDuty(const array<double, 3> d) {
        for (unsigned i = 0; i < 3; ++i) {
            double x = period * (d[i] + 1.0) * 0.5;
            if (x < 0)
                x = 0;
            if (x > period - 1)
                x = period - 1;
            chan[i] = int(x);
        }
    }
    void setSamplePoint(unsigned x) {chan[3] = x;}
    bool operator()(asc::Sampler &sampler, double &t) {
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
static inline void iclarke(const array<T, 2> &ab, array<T, 3> &abc)
{
    abc[0] = ab[0];
    abc[1] = (sqrt(3.0) * ab[1] - ab[0]) / 2.0;
    abc[2] = (-sqrt(3.0) * ab[1] - ab[0]) / 2.0;
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

struct Motor;
class System;
template <typename T>
struct AngleEstimator {
    const T Ts;
    T Vr_a, Vr_b;
    T Va, Vb;
    T Vemf_a, Vemf_b;
    T Vemf_a_ref, Vemf_b_ref;
    T angle;
    T omega;
    T angle_err;
    PIRegS<T> omega_pi;
    PIRegS<T> pi_bemf_a, pi_bemf_b;
    AngleEstimator(const T _Ts) : Ts(_Ts) {
        Vr_a = Vr_b = Vemf_a = Vemf_b = 0;
        Vemf_a_ref = Vemf_b_ref = 0;
        omega = angle = 0.0;
        pi_bemf_a.Ka = pi_bemf_b.Ka = Ls * 4000;
        pi_bemf_a.Kb = pi_bemf_b.Kb = Ts*Rs/Ls;
        omega_pi.Ka = 80;
        omega_pi.Kb = Ts*20;
    }
    void operator()(const array<double, 2> Vab, const array<double, 2> Iab);
};


struct Control {
    array<double, 2> idq;
    array<double, 2> vdq;
    array<double, 2> vab;
    array<double, 3> vabc;
    PIRegulator<double> pi_id;
    PIRegulator<double> pi_iq;
    AngleEstimator<double> a_est;
    bool use_estimator;
    Control(double ts) : idq{0}, vdq{0}, vab{0}, vabc{0}, a_est (ts) {
        pi_iq.Ka = pi_id.Ka = Ls * 2 * M_PI * 4000;
        pi_iq.Kb = pi_id.Kb = ts * Rs/Ls;
        pi_id.limit = pi_iq.limit = VBUS*1.15/2;
        use_estimator = false;
    }
    void operator()(System &sys);
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


/*
static const double Rs = 0.416;
static const double Ls = 1.365e-3;
static const double Ke = 0.2;
static const double poles = 2;
*/
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
    asc::Param Id, Iq;
    asc::Param oe;
    asc::Param theta;
    double Vbus;
    array<double, 3> Vabc;
    array<double, 2> Vab;
    array<double, 2> Vdq;
    array<double, 3> Iabc;
    array<double, 2> Idq;
    array<double, 2> Iab;
    Motor(asc::state_t& state) : Id(state), Iq(state),
        oe(state), theta(state) {
        Iq = Id = 0.0;
        oe = 0.0;
        theta = 0.0;
    }
    void operator()(const asc::state_t&, asc::state_t& D, const double) {
        clarke3(Vabc, Vab);
        park(sin(theta), cos(theta), Vab, Vdq);

        Id(D) = (Vdq[0] + oe * Ls * Iq - Rs * Id) / Ls;
        Iq(D) = (Vdq[1] - oe * (Ls * Id + Ke) - Rs * Iq) / Ls;
        double Te = 1.5 * Ke* poles * Iq;
        double Tf = oe;
        Tf = (Tf > Kf) ?  Kf : Tf;
        Tf = (Tf < -Kf)? -Kf : Tf;
        oe(D) = poles * (Te -Tf - Fv * oe / poles) / (J + Ja);
        theta(D) = oe;
        Idq[0] = Id;
        Idq[1] = Iq;
        ipark(sin(theta), cos(theta), Idq, Iab);
        iclarke(Iab, Iabc);
    }
};

template <typename T>
static inline void svgen(const array<T, 2> &ab, array<T, 3> &abc, const T Vbus)
{
    const T scale = 2.0 / Vbus;
    iclarke(ab, abc);
#if 1
    const T mn = fmin(abc[0], fmin(abc[1], abc[2]));
    const T mx = fmax(abc[0], fmax(abc[1], abc[2]));
    const T off = (mn + mx) * .5;
    for (unsigned i = 0; i < 3; ++i) {
        abc[i] -= off;
        abc[i] *= scale;
    }
#endif
}

struct Simulator;
class System {
    CAPWM pwm;
    Motor motor;
    friend Simulator;
public:
    asc::Param Vca, Vcb, Vcc;
    System(asc::state_t &s) :  pwm(PWM_FREQUENCY, MAX_PERIOD), motor(s), Vca(s), Vcb(s), Vcc(s) {
        motor.Vbus = VBUS;
        Vca = Vcb = Vcc = 0.0;
    }
    double getTheta() {return motor.theta;}
    double getOmega() {return motor.oe;}
    void getCurrents(array<double, 3> &iabc) {iabc = motor.Iabc;}
    double getTs() {return pwm.ts;}
    void operator()(const asc::state_t& x, asc::state_t& D, const double t) {
        for (unsigned i = 0; i < 3; ++i)
            motor.Vabc[i] = motor.Vbus * pwm.Vabc[i];
        motor(x, D, t);

        Vca(D) = ((motor.Vabc[0] - Vca) / R1 - Vca / R2) / C;
        Vcb(D) = ((motor.Vabc[1] - Vcb) / R1 - Vcb / R2) / C;
        Vcc(D) = ((motor.Vabc[2] - Vcc) / R1 - Vcc / R2) / C;
    }
};

template<typename T>
void AngleEstimator<T>::operator()(const array<double, 2> Vab, const array<double, 2> Iab) {
    const double k = (Rs / Ls) * Ts;
    Vr_a += k * (Vab[0] - Vemf_a - Vr_a);
    Vr_b += k * (Vab[1] - Vemf_b - Vr_b);
    Vemf_a = pi_bemf_a(-Iab[0] + Vr_a/Rs);
    Vemf_b = pi_bemf_b(-Iab[1] + Vr_b/Rs);
    angle_err = -2.0*Vemf_a*Vemf_b*cos(2*angle)+
            (Vemf_a*Vemf_a-Vemf_b*Vemf_b)*sin(2*angle);
    omega = omega_pi(angle_err);
    angle += omega * Ts;
    while (angle > 2 * M_PI)
        angle -= 2 * M_PI;
    while (angle < -2 * M_PI)
        angle += 2 * M_PI;
}

void Control::operator()(System &sys)
{
    array<double, 3> Iabc;
    sys.getCurrents(Iabc);
    array<double, 2> Iab, Vab, Idq;
    clarke3(Iabc, Iab);
    clarke3(sys.Vca*VSCALE_DIV, sys.Vcb*VSCALE_DIV, sys.Vcc*VSCALE_DIV, Vab[0], Vab[1]);
    a_est(Vab, Iab);
    double a = sys.getTheta();

    if (use_estimator) {
        a = a_est.angle;
    }
    park(sin(a), cos(a), Iab, Idq);

    // current loop
    vdq[0] = pi_id(idq[0] - Idq[0], -sys.getOmega() * Idq[1] * Ls);
    vdq[1] = pi_iq(idq[1] - Idq[1], sys.getOmega() * (Ls * Idq[0] + Ke));
    ipark(sin(a), cos(a), vdq, vab);
    iclarke(vab, vabc);
    svgen(vab, vabc, VBUS);
}

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
    ctrl(*sys);
    sys->pwm.setDuty(ctrl.vabc);
    ctrl.idq[1] = -2.0;
    while (t < t_end) {
        asc::Sampler sampler(t, dt);
        if (t > 1.0)
            ctrl.idq[1] = 2.0;
        if (t > 1e-2)
            ctrl.use_estimator = true;
        bool tick = false;
        if (sys->pwm(sampler, t)) {
            ctrl(*sys);
            sys->pwm.setDuty(ctrl.vabc);
            tick = true;
        }
        if (tick) rec({t, sys->getOmega(), sys->getTheta()});
        integrator(*sys, s, t, dt);
        while (sys->motor.theta > 2 * M_PI) {
            sys->motor.theta -= 2 * M_PI;
        }
        while (sys->motor.theta < -2 * M_PI) {
            sys->motor.theta += 2 * M_PI;
        }
    }
    rec.csv("pwm", {"t", "omega", "angle"});
}

int main()
{
    Simulator sim;
    Control c(sim.sys->getTs());
    sim.run(c, 2.0);
    return 0;
}
