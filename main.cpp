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

template <typename T>
struct AngleEstimator {
    const T Ts;
    T Vr_a, Vr_b;
    T Vemf_a, Vemf_b;
    T angle;
    T omega, omega_sum;
    T angle_err;
    T K1, K2;
    PIReg<T> pi_bemf_a, pi_bemf_b;
    AngleEstimator(const T _Ts) : Ts(_Ts) {
        Vr_a = Vr_b = Vemf_a = Vemf_b = 0;
        omega_sum = omega = angle = 0.0;
        pi_bemf_a.Kp = pi_bemf_b.Kp = 10.0;
        pi_bemf_a.Ki = pi_bemf_b.Ki = Ts;
        K1 = 1000.0;
        K2 = 1e6;
    }
    void operator()(Motor &motor, T Va, T Vb, T Ia, T Ib);
};


struct Control {
    CAPWM &pwm;
    array<double, 2> idq;
    array<double, 2> vdq;
    array<double, 2> vab;
    array<double, 3> vabc;
    PIRegulator<double> pi_id;
    PIRegulator<double> pi_iq;
    AngleEstimator<double> a_est;
    bool use_estimator;
    Control(CAPWM &_pwm) : pwm(_pwm), idq{0}, vdq{0}, vab{0}, vabc{0}, a_est (pwm.ts) {
        pi_iq.Ka = pi_id.Ka = Ls * 2 * M_PI * 12000;
        pi_iq.Kb = pi_id.Kb = pwm.ts * Rs/Ls;
        pi_id.limit = pi_iq.limit = VBUS*1.15;
        use_estimator = false;
    }
    void operator()(Motor &m);
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
const double C = 47e-9;

struct Motor {
    asc::Param Id, Iq;
    asc::Param Vca, Vcb, Vcc;
    asc::Param oe;
    asc::Param theta;
    Control &ctrl;
    double Vbus;
    array<double, 3> Vabc;
    array<double, 2> Vab;
    array<double, 2> Vdq;
    array<double, 3> Iabc;
    array<double, 2> Idq;
    array<double, 2> Iab;
    Motor(asc::state_t& state, Control &c) : Id(state), Iq(state),
        Vca(state), Vcb(state), Vcc(state), oe(state), theta(state), ctrl(c) {
        Iq = Id = 0.0;
        Vca = Vcb = Vcc = 0.0;
        oe = 0.0;
        theta = 0.0;
    }
    void operator()(const asc::state_t&, asc::state_t& D, const double) {
        for (unsigned i = 0; i < 3; ++i)
            Vabc[i] = Vbus * ctrl.pwm.Vabc[i];
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

        Vca(D) = ((Vabc[0] - Vca) / R1 - Vca / R2) / C;
        Vcb(D) = ((Vabc[1] - Vcb) / R1 - Vcb / R2) / C;
        Vcc(D) = ((Vabc[2] - Vcc) / R1 - Vcc / R2) / C;
    }
};

template<typename T>
void AngleEstimator<T>::operator()(Motor &motor, T Va, T Vb, T Ia, T Ib) {
    //Vr_a += (Va + Vemf_a - Vr_a) * Rs * Ts / Ls;
    //Vr_b += (Vb + Vemf_b - Vr_b) * Rs * Ts / Ls;
    //Vemf_a = pi_bemf_a(Ia - Vr_a/Rs);
    //Vemf_b = pi_bemf_b(Ib - Vr_b/Rs);
    Vemf_a = -motor.oe*Ke*sin(motor.theta);
    Vemf_b = motor.oe*Ke*cos(motor.theta);
    angle_err = -motor.oe*Ke*(sin(angle)*cos(motor.theta) - cos(angle)*sin(motor.theta));
    omega_sum += angle_err * K2 * Ts;
    omega = omega_sum + K1 * angle_err;
    angle += omega * Ts;
    while (angle > 2 * M_PI) {
        angle -= 2*M_PI;
    }
}

const double VSCALE_DIV = (R1 + R2) / R2;

template <typename T>
static inline void svgen(const array<T, 2> &ab, array<T, 3> &abc, const T Vbus)
{
    const T scale = 1.0 / Vbus;
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


void Control::operator()(Motor &m)
{
    array<double, 2> Vab;
    clarke3(m.Vca*VSCALE_DIV, m.Vca*VSCALE_DIV, m.Vca*VSCALE_DIV, Vab[0], Vab[1]);
    a_est(m, Vab[0], Vab[1], m.Iab[0], m.Iab[1]);
    // current loop
    vdq[0] = pi_id(idq[0] - m.Id, -m.oe * m.Iq * Ls);
    vdq[1] = pi_iq(idq[1] - m.Iq, m.oe * (Ls * m.Id + Ke));
    double a = 0;
    if (use_estimator) {
        a = a_est.angle;
    }
    ipark(sin(a), cos(a), vdq, vab);
    iclarke(vab, vabc);
    svgen(vab, vabc, VBUS);
    pwm.setDuty(vabc);
}

int main()
{
    asc::state_t s;
    s.reserve(100);
    CAPWM pwm(PWM_FREQUENCY, MAX_PERIOD);
    Control ctrl(pwm);
    double t = 0.0, t_end = .4, dt = 1e-4;
    asc::RK4 integrator;
    asc::Recorder r1, r2, r3;
    Motor motor(s, ctrl);
    ctrl(motor);
    motor.Vbus = VBUS;
    ctrl.idq[1] = 2.0;
    while (t < t_end) {
        asc::Sampler sampler(t, dt);
        if (t > 1.0)
            ctrl.idq[1] = -2.0;
        if (t > 1e-2)
            ctrl.use_estimator = true;
        if (pwm(sampler, t))
            ctrl(motor);
        r1({t, motor.oe, motor.theta});
        //r1({t, ctrl.vdq[0], ctrl.vdq[1]});
        //t += dt;
        integrator(motor, s, t, dt);
        while (motor.theta > 2 * M_PI) {
            motor.theta -= 2 * M_PI;
        }
        //r2({t, motor.Iabc[0], motor.Iabc[1], motor.Iabc[2]});
        //r2({t, motor.Id, motor.Iq, motor.oe, motor.theta});
        r2({t, -motor.oe*Ke*sin(motor.theta), motor.oe*Ke*cos(motor.theta)});
        r3({t, ctrl.a_est.Vemf_a, ctrl.a_est.Vemf_a,
            ctrl.a_est.omega, ctrl.a_est.angle, ctrl.a_est.angle_err, ctrl.a_est.angle - motor.theta});
    }
    r1.csv("pwm", {"t", "omega", "angle"});
    r2.csv("motor", {"t", "bemf_a", "bemf_b"});
    r3.csv("observer", {"t", "bemf_a", "bemf_b", "omega", "angle", "angle_err", "a_err"});
    return 0;
}
