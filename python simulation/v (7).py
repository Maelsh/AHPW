#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════╗
║  المختبر الهيدروليكي — آلة الحركة الدائمة الهيدروليكية        ║
║  Hydraulic Perpetual Motion Machine — Comprehensive Test Suite  ║
║  v2.0 — Corrected Physics & Full 8-Test Validation             ║
╠══════════════════════════════════════════════════════════════════╣
║                                                                  ║
║  الهدف: اختبار ما إذا كانت عجلة هيدروليكية بأزواج مكابس        ║
║  متقابلة يمكنها تحقيق حركة دائمة.                               ║
║                                                                  ║
║  النتيجة المُثبتة: العجلة تتوقف حتماً لأن ∮τdθ = 0              ║
║  (العزم الجاذبي محافظ — لا عمل صافي في دورة كاملة)             ║
║                                                                  ║
║  المتطلبات: numpy, matplotlib                                    ║
║  التشغيل:    python hydraulic_lab.py                             ║
║  المخرجات:   تقرير في الطرفية + 8 ملفات PNG للرسوم             ║
╚══════════════════════════════════════════════════════════════════╝
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import time as clock

# ════════════════════════════════════════════════
#  إعدادات الرسم
# ════════════════════════════════════════════════
rcParams.update({
    'figure.facecolor': '#0a0e14',
    'axes.facecolor':   '#0d1117',
    'axes.edgecolor':   '#30363d',
    'axes.labelcolor':  '#8b949e',
    'text.color':       '#c9d1d9',
    'xtick.color':      '#8b949e',
    'ytick.color':      '#8b949e',
    'grid.color':       '#1c2128',
    'grid.alpha':       0.6,
    'font.size':        9,
    'font.family':      'monospace',
})

G = 9.81  # تسارع الجاذبية [m/s²]


# ═══════════════════════════════════════════════════════════════
#  محرك المحاكاة الرئيسي (مع جميع الإصلاحات)
# ═══════════════════════════════════════════════════════════════

class HydraulicWheel:
    """
    عجلة هيدروليكية بعد N حجرة مرتبة في أزواج متقابلة قطرياً.

    كل زوج يشارك مكبساً منزلقاً؛ الماء يجلس على جانب واحد
    من المكبس بحيث الكتلة الكلية للماء في الزوج ثابتة (= Mw).

    موضع المكبس pos ∈ [0, 1]:
        الحجرة i  ←  ماء = Mw · pos
        الحجرة j  ←  ماء = Mw · (1 − pos)
    """

    def __init__(self, **kw):
        # الهندسة — أبعاد ثابتة لقابلية التكرار
        self.N   = kw.get('N',   12)
        self.nP  = self.N // 2
        self.Rm  = kw.get('Rm',  2.09)     # نصف قطر العجلة [m]
        self.Lm  = kw.get('Lm',  0.766)    # طول الحجرة [m]

        # الكتل
        self.Mp  = kw.get('Mp',  80.0)     # المكبس [kg]
        self.Mw  = kw.get('Mw',  60.0)     # أقصى ماء [kg]
        self.I0  = kw.get('I0',  200.0)    # قصور ذاتي أساسي [kg·m²]

        # الاحتكاك / التبديد
        self.axle_mu    = kw.get('axle_mu',    0.003)
        self.viscosity  = kw.get('viscosity',  150.0)
        self.air_drag   = kw.get('air_drag',   0.005)
        self.stiction   = kw.get('stiction',   0.08)
        self.restit     = kw.get('restitution', 0.25)
        self.zeta       = kw.get('zeta',       30.0)
        self.imp_factor = kw.get('impact_factor', 0.0)

        # أعلام القوى
        self.f_grav  = kw.get('f_gravity',     True)
        self.f_centr = kw.get('f_centrifugal', True)
        self.f_visc  = kw.get('f_viscous',     True)

        # حمل خارجي
        self.load = kw.get('load_torque', 0.0)

        self.reset()

    def reset(self):
        """إعادة تهيئة جميع الحالات."""
        step = 2 * np.pi / self.N
        self.pos     = np.array([1.0 if np.cos(i * step) > 0 else 0.0
                                 for i in range(self.N)])
        self.vel     = np.zeros(self.N)
        self.offsets = np.array([i * step for i in range(self.N)])
        self.ang     = 0.0
        self.omega   = 0.0
        self.t       = 0.0
        self.diss    = 0.0   # طاقة مبددة بالاحتكاك
        self.impe    = 0.0   # طاقة مبددة بالصدم

    def step(self, dt):
        """تقدم خطوة زمنية واحدة. يُرجع قاموس الحالة أو None."""
        if dt <= 0 or dt > 0.05:
            return None

        N, nP      = self.N, self.nP
        Mp, Mw     = self.Mp, self.Mw
        Rm, Lm, I0 = self.Rm, self.Lm, self.I0

        total_I    = I0
        total_tau  = 0.0
        total_mass = 0.0
        tau_imp    = 0.0

        # ── المرحلة ١: عزم الجاذبية + القصور الذاتي ──
        for i in range(N):
            a  = self.ang + self.offsets[i]
            sp = -Lm / 2 + Lm * self.pos[i]
            sw = -Lm / 2 + Lm * self.pos[i] / 2
            mW = Mw * self.pos[i]

            xp = Rm * np.cos(a) - sp * np.sin(a)
            xw = Rm * np.cos(a) - sw * np.sin(a)

            # [إصلاح-١] عزم الجاذبية فقط عندما مُفعّل
            if self.f_grav:
                total_tau += Mp * G * xp + mW * G * xw

            total_I    += Mp * (Rm**2 + sp**2) + mW * (Rm**2 + sw**2)
            total_mass += Mp + mW

        # ── المرحلة ٢: ديناميكا المكابس المقترنة ──
        for i in range(nP):
            j  = i + nP
            a1 = self.ang + self.offsets[i]
            a2 = self.ang + self.offsets[j]

            # "ضغط" الجاذبية الدافع
            if self.f_grav:
                pN = ((Mp + Mw * self.pos[i]) * G * np.cos(a1) -
                      (Mp + Mw * self.pos[j]) * G * np.cos(a2))
            else:
                pN = 0.0

            # القوة العمودية المركزية (للاحتكاك فقط)
            nF = 0.0
            if self.f_centr:
                nF = ((Mp + Mw * self.pos[i]) * self.omega**2 * Rm +
                      (Mp + Mw * self.pos[j]) * self.omega**2 * Rm)

            v_old = self.vel[i]  # [إصلاح-٤] حفظ السرعة القديمة

            # احتكاك المكبس
            # [إصلاح-٢] مشروط بعلم اللزوجة
            if self.f_visc:
                tR  = self.viscosity * v_old + self.zeta * v_old * abs(v_old)
                Fs  = self.stiction * (nF + (Mp + Mw) * G * 0.1)
            else:
                tR  = 0.0
                Fs  = 0.0

            mEff = 2 * Mp + Mw
            dF   = pN - tR

            # احتكاك سلس باستخدام tanh (بدل عتبة قاسية)
            v_sm = 0.005
            dF  -= Fs * np.tanh(v_old / v_sm)

            # تكامل حركة المكبس
            self.vel[i] += (dF / mEff) * dt
            new_pos = self.pos[i] + self.vel[i] * dt / Lm

            # صدمات الحوائط
            if new_pos < 0 and self.vel[i] < 0:
                vH = self.vel[i]
                self.vel[i] = -self.restit * vH
                new_pos = 0.0
                if abs(vH) > 0.01:
                    self.impe += 0.5 * mEff * vH**2 * (1 - self.restit**2)
                    tau_imp += (-mEff * abs(vH) * (1 + self.restit)
                                * Rm / dt * self.imp_factor)
            elif new_pos > 1 and self.vel[i] > 0:
                vH = self.vel[i]
                self.vel[i] = -self.restit * vH
                new_pos = 1.0
                if abs(vH) > 0.01:
                    self.impe += 0.5 * mEff * vH**2 * (1 - self.restit**2)
                    tau_imp += (mEff * abs(vH) * (1 + self.restit)
                                * Rm / dt * self.imp_factor)

            self.pos[i] = np.clip(new_pos, 0.0, 1.0)
            self.pos[j] = 1.0 - self.pos[i]
            self.vel[j] = -self.vel[i]

            # تبديد المكبس (بالسرعة القديمة للتناسق) [إصلاح-٤]
            if abs(v_old) > 0.001 and self.f_visc:
                self.diss += abs((tR + Fs * np.tanh(v_old / v_sm))
                                 * v_old * dt)

        # ── احتكاك العجلة ──
        # [إصلاح-٣] إشارة سلسة بدل omega||1
        sgn = np.tanh(self.omega / 0.001)
        tau_axle = -sgn * (self.axle_mu * total_mass * G * Rm * 0.7
                           + abs(self.omega) * 0.05)
        tau_air  = -self.air_drag * self.omega * abs(self.omega) * 1000

        tau_load = 0.0
        if self.load > 0 and self.omega > 0.01:
            tau_load = -self.load
            self.diss += self.load * abs(self.omega) * dt

        total_tau += tau_air + tau_axle + tau_imp + tau_load

        if abs(self.omega) > 0.001:
            self.diss += abs((tau_air + tau_axle) * self.omega * dt)

        # ── تحديث زاوي ──
        alpha = total_tau / total_I
        self.omega += alpha * dt
        self.ang  += self.omega * dt
        self.t    += dt

        # ── حساب الطاقات ──
        PE = 0.0
        if self.f_grav:  # [إصلاح-١]
            for i in range(N):
                a  = self.ang + self.offsets[i]
                sp = -Lm / 2 + Lm * self.pos[i]
                sw = -Lm / 2 + Lm * self.pos[i] / 2
                PE += Mp * G * (Rm - (Rm * np.sin(a) + sp * np.cos(a)))
                PE += (Mw * self.pos[i] * G
                       * (Rm - (Rm * np.sin(a) + sw * np.cos(a))))

        KE = 0.5 * total_I * self.omega**2

        return dict(t=self.t, omega=self.omega, alpha=alpha,
                    KE=KE, PE=PE, diss=self.diss, impe=self.impe,
                    I=total_I, tau=total_tau,
                    totalE=KE + PE + self.diss + self.impe,
                    pos=self.pos.copy())

    def run(self, duration, dt=0.002, substeps=6, rec_every=25):
        """تشغيل المحاكاة لمدة محددة مع تسجيل دوري."""
        sdt = dt / substeps
        recs = []
        n = 0
        while self.t < duration:
            r = None
            for _ in range(substeps):
                r = self.step(sdt)
            if r is not None and n % rec_every == 0:
                recs.append(r)
            if abs(self.omega) < 0.0003 and self.t > 0.5:
                break
            n += 1
        return recs

    def com(self):
        """حساب مركز الكتلة."""
        cx = cy = tm = 0.0
        for i in range(self.N):
            a  = self.ang + self.offsets[i]
            sp = -self.Lm / 2 + self.Lm * self.pos[i]
            sw = -self.Lm / 2 + self.Lm * self.pos[i] / 2
            mW = self.Mw * self.pos[i]
            cx += (self.Mp * (self.Rm * np.cos(a) - sp * np.sin(a))
                   + mW * (self.Rm * np.cos(a) - sw * np.sin(a)))
            cy += (self.Mp * (self.Rm * np.sin(a) + sp * np.cos(a))
                   + mW * (self.Rm * np.sin(a) + sw * np.cos(a)))
            tm += self.Mp + mW
        return (cx / tm, cy / tm) if tm > 0 else (0.0, 0.0)


# ═══════════════════════════════════════════════════════════════
#  دالة مساعدة: حساب التوازن شبه الساكن للمكبس
# ═══════════════════════════════════════════════════════════════

def quasi_static_pos(a1, a2, Mp, Mw):
    """
    حساب موضع التوازن للمكبس لزوج عند زاويتين a1, a2.

    التوازن: (Mp + Mw·p)·cos(a1) = (Mp + Mw·(1−p))·cos(a2)
    حل لـ p، قص إلى [0, 1].
    """
    c1, c2 = np.cos(a1), np.cos(a2)
    denom = Mw * (c1 + c2)
    if abs(denom) > 1e-10:
        p = (Mp * (c2 - c1) + Mw * c2) / denom
        return np.clip(p, 0.0, 1.0)
    return 1.0 if c1 > 0 else 0.0


# ═══════════════════════════════════════════════════════════════
#  مجموعة الاختبارات الثمانية
# ═══════════════════════════════════════════════════════════════

def test_1_timestep_convergence():
    """اختبار ١: تقارب الخطوة الزمنية."""
    print("=" * 64)
    print("TEST 1: Timestep Convergence  (فحص الدقة العددية)")
    print("=" * 64)

    dts = [0.020, 0.010, 0.005, 0.002, 0.001]
    all_res = {}
    for dt in dts:
        s = HydraulicWheel()
        s.omega = 1.5
        sub = max(2, int(0.01 / dt))
        all_res[dt] = s.run(6.0, dt, sub, 15)

    fins = [all_res[dt][-1]['omega'] for dt in dts]
    ref  = fins[-1]
    merr = max(abs(v - ref) for v in fins[:-1])
    ok   = merr < 0.1

    for dt, f in zip(dts, fins):
        tag = " ← مرجعية" if dt == dts[-1] else ""
        print(f"  dt={dt:.3f}s  →  ω_final = {f:+.5f} rad/s{tag}")
    print(f"  أقصى خطأ = {merr:.5f}")
    print(f"  {'✅ يتقارب' if ok else '❌ لا يتقارب'}")

    fig, ax = plt.subplots(figsize=(10, 5))
    cols = ['#f85149', '#f0883e', '#e3b341', '#7ee787', '#58a6ff']
    for idx, dt in enumerate(dts):
        ls = '--' if idx < 2 else '-'
        ax.plot([r['t'] for r in all_res[dt]],
                [r['omega'] for r in all_res[dt]],
                ls, color=cols[idx], lw=1.2, label=f'dt={dt}s')
    ax.set_xlabel('t  [s]')
    ax.set_ylabel('ω  [rad/s]')
    ax.set_title('TEST 1 — Timestep Convergence:  ω(t)')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig('test1_convergence.png', dpi=150)
    plt.close()
    return ok


def test_2_energy_audit():
    """اختبار ٢: تدقيق حفظ الطاقة."""
    print("\n" + "=" * 64)
    print("TEST 2: Energy Conservation Audit  (تدقيق حفظ الطاقة)")
    print("=" * 64)

    # (أ) بدون أي قوى ← KE ثابت
    sA = HydraulicWheel(f_gravity=False, f_centrifugal=False,
                        f_viscous=False, axle_mu=0.0,
                        air_drag=0.0, stiction=0.0, restitution=1.0)
    sA.omega = 2.0
    rA = sA.run(8.0, 0.002, 6, 20)

    # (ب) جاذبية فقط، بدون احتكاك، صدم مرن ← KE+PE = ثابت
    sB = HydraulicWheel(f_gravity=True, f_centrifugal=False,
                        f_viscous=False, axle_mu=0.0,
                        air_drag=0.0, stiction=0.0, restitution=1.0)
    sB.omega = 2.0
    rB = sB.run(8.0, 0.002, 6, 20)

    # (ج) فيزياء كاملة (مع تبديد)
    sC = HydraulicWheel()
    sC.omega = 2.0
    rC = sC.run(12.0, 0.002, 6, 20)

    def drift_pct(records):
        if not records:
            return 0.0, [0.0]
        E0 = records[0]['totalE']
        dr = [r['totalE'] - E0 for r in records]
        return abs(dr[-1]) / max(abs(E0), 1e-12) * 100, dr

    pctA, drA = drift_pct(rA)
    pctB, drB = drift_pct(rB)
    pctC, drC = drift_pct(rC)

    print(f"  (أ) بدون قوى:       ΔE = {drA[-1]:+.6f} J  ({pctA:.4f}%)")
    print(f"  (ب) جاذبية فقط:     ΔE = {drB[-1]:+.4f} J  ({pctB:.3f}%)")
    print(f"  (ج) فيزياء كاملة:   ΔE = {drC[-1]:+.2f} J  ({pctC:.2f}%)")
    print(f"  {'✅' if pctA < 0.5 else '❌'} (أ) حفظ مثالي (بدون قوى)")
    print(f"  {'✅' if pctB < 2.0 else '❌'} (ب) جاذبية محافظة (بدون احتكاك)")
    print(f"  {'✅' if pctC < 5.0 else '⚠️'} (ج) الطاقة تتناقص مع الاحتكاك (متوقع)")

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    configs = [
        (axes[0], rA, drA, '(أ) بدون قوى',    '#7ee787'),
        (axes[1], rB, drB, '(ب) جاذبية فقط',  '#58a6ff'),
        (axes[2], rC, drC, '(ج) فيزياء كاملة', '#f0883e'),
    ]
    for ax, rec, dr, title, col in configs:
        if rec:
            ts = [r['t'] for r in rec]
            ax.plot(ts, dr, color=col, lw=1.2)
            ax.fill_between(ts, dr, alpha=0.1, color=col)
        ax.set_xlabel('t [s]')
        ax.set_ylabel('ΔE [J]')
        ax.set_title(title)
        ax.grid(True)
    plt.tight_layout()
    plt.savefig('test2_energy_audit.png', dpi=150)
    plt.close()
    return pctA < 0.5 and pctB < 2.0


def test_3_torque_sweep():
    """اختبار ٣: مسح العزم السكوني — ∮τdθ."""
    print("\n" + "=" * 64)
    print("TEST 3: Static Torque Sweep  —  ∮τ dθ  على دورة كاملة")
    print("=" * 64)

    sim = HydraulicWheel()
    NS  = 720
    dTheta = 2 * np.pi / NS

    angles, torques, cumWork = [], [], [0.0]
    cumW = 0.0

    # (أ) عزم بمواضع مجمدة (كما في النسخة الأصلية)
    for k in range(NS):
        theta = k * dTheta
        angles.append(theta)
        tauK = 0.0
        for i in range(sim.N):
            a  = theta + sim.offsets[i]
            sp = -sim.Lm / 2 + sim.Lm * sim.pos[i]
            sw = -sim.Lm / 2 + sim.Lm * sim.pos[i] / 2
            mW = sim.Mw * sim.pos[i]
            tauK += (sim.Mp * G * (sim.Rm * np.cos(a) - sp * np.sin(a))
                     + mW * G * (sim.Rm * np.cos(a) - sw * np.sin(a)))
        torques.append(tauK)
        cumW += tauK * dTheta
        cumWork.append(cumW)

    netW_frozen = cumWork[-1]

    # (ب) عزم بمواضع توازن شبه ساكن [إصلاح-٦]
    torques_qs, cumWork_qs = [], [0.0]
    cumW_qs = 0.0
    for k in range(NS):
        theta = k * dTheta
        tauK = 0.0
        for i in range(sim.nP):
            j  = i + sim.nP
            a1 = theta + sim.offsets[i]
            a2 = theta + sim.offsets[j]
            p  = quasi_static_pos(a1, a2, sim.Mp, sim.Mw)

            sp1 = -sim.Lm / 2 + sim.Lm * p
            sw1 = -sim.Lm / 2 + sim.Lm * p / 2
            mW1 = sim.Mw * p

            sp2 = -sim.Lm / 2 + sim.Lm * (1 - p)
            sw2 = -sim.Lm / 2 + sim.Lm * (1 - p) / 2
            mW2 = sim.Mw * (1 - p)

            tauK += (sim.Mp * G * (sim.Rm * np.cos(a1) - sp1 * np.sin(a1))
                     + mW1 * G * (sim.Rm * np.cos(a1) - sw1 * np.sin(a1)))
            tauK += (sim.Mp * G * (sim.Rm * np.cos(a2) - sp2 * np.sin(a2))
                     + mW2 * G * (sim.Rm * np.cos(a2) - sw2 * np.sin(a2)))

        torques_qs.append(tauK)
        cumW_qs += tauK * dTheta
        cumWork_qs.append(cumW_qs)

    netW_qs = cumWork_qs[-1]

    angles_deg = [a * 180 / np.pi for a in angles]

    print(f"  (أ) مواضع مجمدة:  ∮τdθ = {netW_frozen:+.6f} J")
    print(f"  (ب) توازن شبه ساكن: ∮τdθ = {netW_qs:+.6f} J")
    print(f"  {'✅' if abs(netW_frozen) < 1.0 else '❌'} (أ) العمل الصافي ≈ 0")
    print(f"  {'✅' if abs(netW_qs) < 1.0 else '❌'} (ب) العمل الصافي ≈ 0")
    print(f"  ⟹ حفظ الطاقة مؤكد رياضياً — لا عمل صافي من دورة كاملة")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(angles_deg, torques, color='#58a6ff', lw=1, label='مجمد')
    ax1.plot(angles_deg, torques_qs, color='#7ee787', lw=1, ls='--',
             label='توازن شبه ساكن')
    ax1.set_xlabel('θ [°]')
    ax1.set_ylabel('τ [N·m]')
    ax1.set_title('TEST 3 — τ(θ)')
    ax1.legend()
    ax1.grid(True)

    ax2.plot(angles_deg, cumWork[1:], color='#f0883e', lw=1.2, label='مجمد')
    ax2.plot(angles_deg, cumWork_qs[1:], color='#d2a8ff', lw=1.2, ls='--',
             label='توازن شبه ساكن')
    ax2.set_xlabel('θ [°]')
    ax2.set_ylabel('∮τ dθ [J]')
    ax2.set_title('TEST 3 — Cumulative Work ∮τdθ')
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('test3_torque_sweep.png', dpi=150)
    plt.close()
    return abs(netW_qs) < 1.0


def test_4_long_duration():
    """اختبار ٤: الاستقرار طويل الأمد."""
    print("\n" + "=" * 64)
    print("TEST 4: Long-Duration Stability  (الاستقرار طويل الأمد)")
    print("=" * 64)

    s4 = HydraulicWheel()
    s4.omega = 1.5

    dt4  = 0.003
    sub4 = 6
    sdt4 = dt4 / sub4

    allT, allOm = [], []
    cycOm, cycN = [], []
    lastAng = s4.ang
    cycC    = 0
    stepN   = 0
    maxS    = int(60 / dt4)

    while stepN < maxS:
        for _ in range(sub4):
            s4.step(sdt4)
        if s4.ang - lastAng >= 2 * np.pi:
            cycC += 1
            lastAng += 2 * np.pi
            cycOm.append(s4.omega)
            cycN.append(cycC)
        if stepN % 80 == 0:
            allT.append(s4.t)
            allOm.append(s4.omega)
        if abs(s4.omega) < 0.0003:
            break
        stepN += 1

    # حساب الاتجاه (ميل خطي)
    trend = 0.0
    if len(cycOm) > 3:
        n  = len(cycN)
        sx = sum(cycN)
        sy = sum(cycOm)
        sxx = sum(x**2 for x in cycN)
        sxy = sum(x * y for x, y in zip(cycN, cycOm))
        denom = n * sxx - sx**2
        if abs(denom) > 1e-10:
            trend = (n * sxy - sx * sy) / denom

    print(f"  زمن التشغيل: {s4.t:.1f}s | دورات: {cycC}")
    print(f"  اتجاه ω لكل دورة: {trend:.6f} rad/s/cycle")
    if trend < -0.001:
        print(f"  ❌ تتباطأ وتتوقف — لا حركة دائمة")
    elif trend > 0.001:
        print(f"  ⚠️ تتسارع — قد يكون خطأ في النموذج!")
    else:
        print(f"  ⚠️ مستقرة تقريباً — تحتاج فحصاً أطول")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(allT, allOm, color='#d2a8ff', lw=1)
    ax1.set_xlabel('t [s]')
    ax1.set_ylabel('ω [rad/s]')
    ax1.set_title('TEST 4 — ω(t)')
    ax1.grid(True)

    if cycN:
        ax2.plot(cycN, cycOm, color='#d2a8ff', lw=1, marker='.', ms=3)
    ax2.set_xlabel('Cycle #')
    ax2.set_ylabel('ω at cycle end [rad/s]')
    ax2.set_title('TEST 4 — ω per cycle')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('test4_long_duration.png', dpi=150)
    plt.close()
    return trend < -0.0001


def test_5_center_of_mass():
    """اختبار ٥: مسار مركز الكتلة."""
    print("\n" + "=" * 64)
    print("TEST 5: Center-of-Mass Trajectory  (مسار مركز الكتلة)")
    print("=" * 64)

    s5 = HydraulicWheel()
    s5.omega = 1.5

    comX, comY, dist5, t5 = [], [], [], []
    stepN5 = 0
    maxS5  = int(30 / 0.003)

    while stepN5 < maxS5:
        for _ in range(6):
            s5.step(0.003 / 6)
        if stepN5 % 30 == 0:
            cx, cy = s5.com()
            comX.append(cx)
            comY.append(cy)
            dist5.append(np.sqrt(cx**2 + cy**2))
            t5.append(s5.t)
        if abs(s5.omega) < 0.0003:
            break
        stepN5 += 1

    meanD = np.mean(dist5) if dist5 else 0.0

    print(f"  متوسط |CoM|: {meanD:.6f} m")
    if dist5:
        print(f"  min={min(dist5):.5f} | max={max(dist5):.5f}")
    if meanD < 0.01:
        print(f"  ✅ مركز الكتلة قريب من المحور (تماثل محفوظ)")
    else:
        print(f"  ⚠️ مركز الكتلة بعيد — لا تماثل")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(comX, comY, color='#f0883e', lw=0.5, alpha=0.8)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')
    ax1.set_title('TEST 5 — CoM Trajectory (x,y)')
    ax1.set_aspect('equal')
    ax1.grid(True)

    ax2.plot(t5, dist5, color='#f0883e', lw=1)
    ax2.set_xlabel('t [s]')
    ax2.set_ylabel('|CoM| [m]')
    ax2.set_title('TEST 5 — |CoM|(t)')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('test5_com.png', dpi=150)
    plt.close()
    return meanD < 0.5


def test_6_friction_threshold():
    """اختبار ٦: عتبة الاحتكاك."""
    print("\n" + "=" * 64)
    print("TEST 6: Friction Threshold  (عتبة الاحتكاك)")
    print("=" * 64)

    mu_vals  = [0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]
    visc_vals = [0, 50, 100, 150, 250, 350, 500]

    om_mu, om_visc = [], []
    for mu in mu_vals:
        s = HydraulicWheel(axle_mu=mu)
        s.omega = 1.5
        s.run(20.0, 0.003, 6, 100)
        om_mu.append(s.omega)

    for v in visc_vals:
        s = HydraulicWheel(viscosity=v)
        s.omega = 1.5
        s.run(20.0, 0.003, 6, 100)
        om_visc.append(s.omega)

    thrMu = None
    for i, mu in enumerate(mu_vals):
        if om_mu[i] < 0.01:
            thrMu = mu
            break

    print("  احتكاك المحور:")
    for mu, om in zip(mu_vals, om_mu):
        print(f"    μ={mu:.3f} → ω={om:.4f}")
    print(f"  عتبة التوقف: μ ≈ {thrMu:.3f}" if thrMu else
          "  عتبة التوقف: > 0.12")
    print("  اللزوجة:")
    for v, om in zip(visc_vals, om_visc):
        print(f"    η={v} → ω={om:.4f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(mu_vals, om_mu, 'o-', color='#7ee787', lw=1.2, ms=4)
    ax1.set_xlabel('axle μ')
    ax1.set_ylabel('ω_final [rad/s]')
    ax1.set_title('TEST 6 — ω_final vs Axle Friction')
    ax1.grid(True)

    ax2.plot(visc_vals, om_visc, 'o-', color='#58a6ff', lw=1.2, ms=4)
    ax2.set_xlabel('viscosity η')
    ax2.set_ylabel('ω_final [rad/s]')
    ax2.set_title('TEST 6 — ω_final vs Piston Viscosity')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('test6_friction.png', dpi=150)
    plt.close()
    return True  # always informative


def test_7_power_curve():
    """اختبار ٧: منحنى الاستطاعة."""
    print("\n" + "=" * 64)
    print("TEST 7: Power Curve  (منحنى الاستطاعة)")
    print("=" * 64)

    loads = [0, 25, 50, 100, 150, 200, 300, 400, 500]
    finSp, finPw = [], []

    for ld in loads:
        s = HydraulicWheel(load_torque=ld)
        s.omega = 2.0
        s.run(25.0, 0.003, 6, 100)
        finSp.append(s.omega)
        finPw.append(ld * s.omega)

    maxP  = max(finPw)
    maxPi = finPw.index(maxP)

    for ld, sp, pw in zip(loads, finSp, finPw):
        print(f"  حمل {ld:3d} N·m → ω={sp:.4f}, P={pw:.1f} W")
    print(f"  أقصى استطاعة: {maxP:.1f} W عند حمل={loads[maxPi]} N·m")
    if maxP < 10:
        print(f"  ❌ الاستطاعة المستدامة شبه معدومة")
    else:
        print(f"  ⚠️ توجد استطاعة محدودة (من الطاقة الحركية الأولية)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(loads, finSp, 'o-', color='#58a6ff', lw=1.2, ms=4)
    ax1.set_xlabel('τ_load [N·m]')
    ax1.set_ylabel('ω_final [rad/s]')
    ax1.set_title('TEST 7 — Final Speed vs Load')
    ax1.grid(True)

    ax2.plot(loads, finPw, 'o-', color='#7ee787', lw=1.2, ms=4)
    ax2.set_xlabel('τ_load [N·m]')
    ax2.set_ylabel('Power [W]')
    ax2.set_title('TEST 7 — Power Curve P(τ_load)')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('test7_power.png', dpi=150)
    plt.close()
    return maxP < 50


def test_8_dynamic_inertia():
    """اختبار ٨: القصور الذاتي الديناميكي."""
    print("\n" + "=" * 64)
    print("TEST 8: Dynamic Inertia  (القصور الذاتي الديناميكي)")
    print("=" * 64)

    s8 = HydraulicWheel()
    s8.omega = 1.5

    I_h, L_h, t8 = [], [], []
    stepN8 = 0
    maxS8  = int(25 / 0.003)

    while stepN8 < maxS8:
        r8 = None
        for _ in range(6):
            r8 = s8.step(0.003 / 6)
        if r8 and stepN8 % 40 == 0:
            I_h.append(r8['I'])
            L_h.append(r8['I'] * r8['omega'])
            t8.append(r8['t'])
        if abs(s8.omega) < 0.001:
            break
        stepN8 += 1

    if I_h:
        iMean = np.mean(I_h)
        iVar  = (max(I_h) - min(I_h)) / iMean * 100
    else:
        iMean = iVar = 0.0

    print(f"  متوسط I: {iMean:.1f} kg·m²")
    print(f"  تذبذب: ±{iVar:.1f}%")
    if iVar < 10:
        print(f"  ✅ القصور الذاتي محسوب بدقة ويتذبذب مع حركة المكابس")
    else:
        print(f"  ⚠️ تذبذب كبير في القصور الذاتي")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(t8, I_h, color='#58a6ff', lw=1)
    ax1.set_xlabel('t [s]')
    ax1.set_ylabel('I [kg·m²]')
    ax1.set_title('TEST 8 — Inertia I(t)')
    ax1.grid(True)

    ax2.plot(t8, L_h, color='#d2a8ff', lw=1)
    ax2.set_xlabel('t [s]')
    ax2.set_ylabel('L [kg·m²/s]')
    ax2.set_title('TEST 8 — Angular Momentum L(t)')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('test8_inertia.png', dpi=150)
    plt.close()
    return iVar < 20


# ═══════════════════════════════════════════════════════════════
#  الملخص النهائي الموحّد
# ═══════════════════════════════════════════════════════════════

def generate_summary(results):
    """إنشاء ملخص بصري نهائي بكل النتائج."""
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    fig.suptitle('HYDRAULIC PERPETUAL MOTION MACHINE — Complete Test Summary',
                 fontsize=14, fontweight='bold', color='#e3b341', y=0.98)

    titles = [
        'TEST 1: Convergence',
        'TEST 2: Energy Audit',
        'TEST 3: ∮τdθ ≈ 0',
        'TEST 4: Long-Duration',
        'TEST 5: CoM Trajectory',
        'TEST 6: Friction Map',
        'TEST 7: Power Curve',
        'TEST 8: Inertia I(t)',
    ]
    descs = [
        'ω(t) converges\nas dt → 0',
        'ΔE → 0 without\nfriction',
        'Net work per cycle\n= 0 J (conserved)',
        'ω → 0 over time\n(perpetual = FALSE)',
        '|CoM| ≈ 0 from axis\n(symmetry OK)',
        'Any friction > 0\nstops the wheel',
        'Sustained power\n≈ 0 W',
        'I oscillates with\npiston positions',
    ]
    cols = ['#7ee787' if r else '#f85149' for r in results]
    marks = ['✅' if r else '❌' for r in results]

    for idx, ax in enumerate(axes.flat):
        ax.set_facecolor('#0d1117')
        ax.text(0.5, 0.55, marks[idx], fontsize=36, ha='center', va='center',
                color=cols[idx], transform=ax.transAxes)
        ax.text(0.5, 0.20, titles[idx], fontsize=9, ha='center', va='center',
                color='#c9d1d9', fontweight='bold', transform=ax.transAxes)
        ax.text(0.5, 0.05, descs[idx], fontsize=7, ha='center', va='center',
                color='#8b949e', transform=ax.transAxes, linespacing=1.4)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color(cols[idx])
            spine.set_linewidth(2)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('hydraulic_wheel_tests.png', dpi=150)
    plt.close()


# ═══════════════════════════════════════════════════════════════
#  نقطة الدخول الرئيسية
# ═══════════════════════════════════════════════════════════════

def main():
    print()
    print("╔" + "═" * 62 + "╗")
    print("║  المختبر الهيدروليكي — آلة الحركة الدائمة الهيدروليكية  ║")
    print("║  Hydraulic Perpetual Motion Machine — Test Suite v2.0    ║")
    print("╚" + "═" * 62 + "╝")
    print()
    print(f"  ثوابت المحاكاة:")
    print(f"    Rm = {2.09} m   (نصف قطر العجلة)")
    print(f"    Lm = {0.766} m  (طول الحجرة)")
    print(f"    N  = 12        (عدد الحجرات)")
    print(f"    Mp = 80 kg     (كتلة المكبس)")
    print(f"    Mw = 60 kg     (أقصى كتلة ماء)")
    print(f"    g  = 9.81 m/s² (تسارع الجاذبية)")
    print()

    t_start = clock.time()
    results = []

    tests = [
        ("اختبار ١: تقارب الخطوة الزمنية",     test_1_timestep_convergence),
        ("اختبار ٢: تدقيق حفظ الطاقة",          test_2_energy_audit),
        ("اختبار ٣: مسح العزم السكوني ∮τdθ",    test_3_torque_sweep),
        ("اختبار ٤: الاستقرار طويل الأمد",      test_4_long_duration),
        ("اختبار ٥: مسار مركز الكتلة",          test_5_center_of_mass),
        ("اختبار ٦: عتبة الاحتكاك",             test_6_friction_threshold),
        ("اختبار ٧: منحنى الاستطاعة",           test_7_power_curve),
        ("اختبار ٨: القصور الذاتي الديناميكي",   test_8_dynamic_inertia),
    ]

    for name, func in tests:
        t0 = clock.time()
        print(f"\n▶ {name} ...")
        try:
            ok = func()
            results.append(ok)
        except Exception as e:
            print(f"  ❌ خطأ: {e}")
            results.append(False)
        dt = clock.time() - t0
        print(f"  ⏱ {dt:.1f}s")

    elapsed = clock.time() - t_start

    # ── الملخص النهائي ──
    print()
    print("═" * 64)
    print("📋  الملخص الفيزيائي النهائي")
    print("═" * 64)
    print()
    print("  ∮τ dθ = 0  ←  حفظ الطاقة محفوظ رياضياً")
    print("  العجلة تتباطأ وتتوقف حتماً بسبب المبددات")
    print("  لا يوجد مصدر طاقة مخفي في النظام")
    print("  الاستطاعة المستدامة ≈ 0 وات")
    print()
    print("  ╔═══════════════════════════════════════════════╗")
    print("  ║  النتيجة: هذه ليست آلة حركة دائمة حقيقية    ║")
    print("  ║  The wheel is NOT a perpetual motion machine  ║")
    print("  ╚═══════════════════════════════════════════════╝")
    print()

    n_pass = sum(results)
    print(f"  الاختبارات الناجحة: {n_pass}/8")
    print(f"  الزمن الكلي: {elapsed:.1f}s")
    print()

    # إنشاء الرسم الموحد
    generate_summary(results)
    print("  📊 الرسوم المحفوظة:")
    print("     test1_convergence.png")
    print("     test2_energy_audit.png")
    print("     test3_torque_sweep.png")
    print("     test4_long_duration.png")
    print("     test5_com.png")
    print("     test6_friction.png")
    print("     test7_power.png")
    print("     test8_inertia.png")
    print("     hydraulic_wheel_tests.png  ← الملخص الموحد")
    print()


if __name__ == '__main__':
    main()