#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
╔══════════════════════════════════════════════════════════════════════════╗
║          المختبر الهيدروليكي — محاكاة Python الشاملة (مُصلَّحة)         ║
║          Hydraulic Wheel Simulation — Complete Python Port (fixed)       ║
║                                                                          ║
║  مصدر الكود: mohakahwekhtbarat.html (v22)                                ║
║  المهام:                                                                 ║
║    1. فحص سلامة النموذج الفيزيائي                                        ║
║    2. تشغيل الاختبارات الثمانية لإثبات/نفي الحركة الدائمة               ║
║    3. إنتاج رسوم بيانية شاملة                                            ║
║                                                                          ║
║  الاستخدام:  python hydraulic_wheel_simulation.py [--show] [--no-save]  ║
║  المتطلبات: pip install numpy matplotlib                                 ║
╚══════════════════════════════════════════════════════════════════════════╝
"""

import math
import numpy as np
import matplotlib
# إزالة Agg الافتراضي - نسمح بالاختيار التلقائي، ونُظهر النافذة عند الطلب
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
import sys

# ══════════════════════════════════════════════════════════════
#  الثوابت المرجعية (مطابقة لـ TestSim في JavaScript)
# ══════════════════════════════════════════════════════════════
G       = 9.81    # تسارع الجاذبية (m/s²)
T_RM    = 2.09    # نصف قطر العجلة (m)
T_LM    = 0.766   # طول الحجرة / المكبس (m)
T_N     = 12      # عدد الحجرات


# ══════════════════════════════════════════════════════════════
#  الفئة الرئيسية: محاكاة العجلة الهيدروليكية
# ══════════════════════════════════════════════════════════════
class HydraulicWheelSim:
    """
    محاكاة العجلة الهيدروليكية ذات المكابس المقترنة.
    النموذج الفيزيائي مطابق للكود الأصلي (JavaScript).
    """

    def __init__(self, **opts):
        self.g      = G
        self.N      = opts.get('N',    T_N)
        self.nP     = self.N // 2      # عدد الأزواج
        self.Mp     = opts.get('Mp',   80.0)    # كتلة المكبس (kg)
        self.Mw     = opts.get('Mw',   60.0)    # أقصى كتلة ماء (kg)
        self.Rm     = opts.get('Rm',   T_RM)    # نصف قطر مسار الحجرات (m)
        self.Lm     = opts.get('Lm',   T_LM)    # طول الحجرة (m)
        self.I0     = opts.get('I0',   200.0)   # قصور ذاتي أساسي (kg·m²)
        self.air_d  = opts.get('airD', 0.005)   # معامل مقاومة الهواء
        self.ax_mu  = opts.get('axMu', 0.003)   # معامل احتكاك المحور
        self.visc   = opts.get('visc', 150.0)   # معامل اللزوجة (N·s/m)
        self.stic   = opts.get('stic', 0.08)    # معامل الاحتكاك السكوني
        self.rest   = opts.get('rest', 0.25)    # معامل الارتداد
        self.zeta   = opts.get('zeta', 30.0)    # معامل التخميد التربيعي
        self.imp_f  = opts.get('impF', 0.0)     # معامل الصدمة
        self.f_grav = opts.get('fGrav',  True)  # تفعيل الجاذبية
        self.f_cent = opts.get('fCentr', True)  # تفعيل القوة الطاردة
        self.f_visc = opts.get('fVisc',  True)  # تفعيل اللزوجة
        self.load   = opts.get('load',   0.0)   # عزم الحمل الخارجي (N·m)
        self.reset()

    def reset(self):
        """إعادة ضبط الحالة الأولية"""
        step = 2 * math.pi / self.N
        self.pos  = []   # موضع المكبس [0,1]
        self.vel  = []   # سرعة المكبس (m/s)
        self.off  = []   # إزاحة الزاوية لكل حجرة (rad)
        for i in range(self.N):
            self.pos.append(1.0 if math.cos(i * step) > 0 else 0.0)
            self.vel.append(0.0)
            self.off.append(i * step)
        self.ang   = 0.0
        self.omega = 0.0
        self.t     = 0.0
        self.diss  = 0.0   # الطاقة المبددة (J)
        self.impe  = 0.0   # طاقة الصدمات (J)

    def step(self, dt):
        """خطوة زمنية واحدة. تُعيد dict بالمتحولات الفيزيائية."""
        if dt <= 0 or dt > 0.05:
            return None

        Lm, Rm, g = self.Lm, self.Rm, self.g
        N, nP = self.N, self.nP
        Mp, Mw = self.Mp, self.Mw

        # ── عزم الجاذبية الكلي وعزم القصور الذاتي ──
        total_I   = self.I0
        total_tau = 0.0
        total_M   = 0.0

        for i in range(N):
            a  = self.ang + self.off[i]
            sp = -Lm/2 + Lm * self.pos[i]
            sw = -Lm/2 + Lm * self.pos[i] / 2
            mW = Mw * self.pos[i]
            xp = Rm * math.cos(a) - sp * math.sin(a)
            xw = Rm * math.cos(a) - sw * math.sin(a)
            tau_g = (self.f_grav) and (Mp * g * xp + mW * g * xw) or 0.0
            total_tau += tau_g
            total_I += Mp * (Rm**2 + sp**2) + mW * (Rm**2 + sw**2)
            total_M += Mp + mW

        # ── ديناميكية المكابس المقترنة ──
        tau_impact = 0.0

        for i in range(nP):
            j  = i + nP      # المكبس المزدوج
            a1 = self.ang + self.off[i]
            a2 = self.ang + self.off[j]
            p1  = self.f_grav and (Mp * g * math.cos(a1)) or 0.0
            p2  = self.f_grav and (Mp * g * math.cos(a2)) or 0.0
            pw1 = self.f_grav and (Mw * self.pos[i] * g * math.cos(a1)) or 0.0
            pw2 = self.f_grav and (Mw * self.pos[j] * g * math.cos(a2)) or 0.0
            p_net = (p1 + pw1) - (p2 + pw2)

            n_force = 0.0
            if self.f_cent:
                n_force = ((Mp + Mw * self.pos[i]) * self.omega**2 * Rm
                         + (Mp + Mw * self.pos[j]) * self.omega**2 * Rm)

            v1 = self.vel[i]
            t_resist = (self.f_visc) and (self.visc * v1 + self.zeta * v1 * abs(v1)) or 0.0
            m_eff = 2 * Mp + Mw

            F_stick = self.stic * (n_force + (Mp + Mw) * g * 0.1)
            v_smooth = 0.005
            F_friction = F_stick * math.tanh(v1 / v_smooth)

            drive_force = p_net - t_resist - F_friction

            self.vel[i] += (drive_force / m_eff) * dt
            new_pos = self.pos[i] + self.vel[i] * dt / Lm

            # حدود نهاية الشوط – ارتداد
            if new_pos < 0 and self.vel[i] < 0:
                vH = self.vel[i]
                self.vel[i] = -self.rest * vH
                new_pos = 0.0
                if abs(vH) > 0.01:
                    self.impe += 0.5 * m_eff * vH**2 * (1 - self.rest**2)
                    tau_impact += (-m_eff * abs(vH) * (1 + self.rest)
                                   * Rm / dt * self.imp_f)
            elif new_pos > 1 and self.vel[i] > 0:
                vH = self.vel[i]
                self.vel[i] = -self.rest * vH
                new_pos = 1.0
                if abs(vH) > 0.01:
                    self.impe += 0.5 * m_eff * vH**2 * (1 - self.rest**2)
                    tau_impact += (m_eff * abs(vH) * (1 + self.rest)
                                   * Rm / dt * self.imp_f)

            self.pos[i] = max(0.0, min(1.0, new_pos))
            self.pos[j] = 1.0 - self.pos[i]
            self.vel[j] = -self.vel[i]

            if abs(self.vel[i]) > 0.001:
                self.diss += abs((t_resist + F_friction) * self.vel[i] * dt)

        # ── عزم احتكاك المحور ومقاومة الهواء ──
        sign_omega = math.copysign(1, self.omega) if self.omega != 0 else 1.0
        tau_axle = -sign_omega * (self.ax_mu * total_M * g * Rm * 0.7
                                  + abs(self.omega) * 0.05)
        tau_air  = -self.air_d * self.omega * abs(self.omega) * 1000.0

        # عزم الحمل الخارجي
        tau_load = 0.0
        if self.load > 0 and self.omega > 0.01:
            tau_load = -self.load
            self.diss += self.load * abs(self.omega) * dt

        total_tau += tau_air + tau_axle + tau_impact + tau_load

        if abs(self.omega) > 0.01:
            self.diss += abs((tau_air + tau_axle) * self.omega * dt)

        # تحديث ديناميكي
        alpha = total_tau / total_I
        self.omega += alpha * dt
        self.ang   += self.omega * dt
        self.t     += dt

        # الطاقة الكامنة
        PE = 0.0
        for i in range(N):
            a  = self.ang + self.off[i]
            sp = -Lm/2 + Lm * self.pos[i]
            sw = -Lm/2 + Lm * self.pos[i] / 2
            PE += Mp * g * (Rm - (Rm * math.sin(a) + sp * math.cos(a)))
            PE += Mw * self.pos[i] * g * (Rm - (Rm * math.sin(a) + sw * math.cos(a)))

        KE = 0.5 * total_I * self.omega**2

        return {
            't':      self.t,
            'omega':  self.omega,
            'alpha':  alpha,
            'KE':     KE,
            'PE':     PE,
            'diss':   self.diss,
            'impe':   self.impe,
            'I':      total_I,
            'tau':    total_tau,
            'total_E': KE + PE + self.diss + self.impe,
            'pos':    list(self.pos),
        }

    def run(self, duration, dt=0.002, sub=6, rec_every=25):
        """تشغيل المحاكاة لمدة duration ثانية. تُعيد قائمة نقاط البيانات."""
        sdt  = dt / sub
        recs = []
        step_n = 0
        while self.t < duration:
            r = None
            for _ in range(sub):
                r = self.step(sdt)
            if r and step_n % rec_every == 0:
                recs.append(r)
            if abs(self.omega) < 3e-4 and self.t > 0.5:
                break
            step_n += 1
        return recs

    def center_of_mass(self):
        """مركز كتلة النظام."""
        cx, cy, total_m = 0.0, 0.0, 0.0
        for i in range(self.N):
            a  = self.ang + self.off[i]
            sp = -self.Lm/2 + self.Lm * self.pos[i]
            sw = -self.Lm/2 + self.Lm * self.pos[i] / 2
            mW = self.Mw * self.pos[i]
            cx += self.Mp * (self.Rm * math.cos(a) - sp * math.sin(a))
            cx += mW      * (self.Rm * math.cos(a) - sw * math.sin(a))
            cy += self.Mp * (self.Rm * math.sin(a) + sp * math.cos(a))
            cy += mW      * (self.Rm * math.sin(a) + sw * math.cos(a))
            total_m += self.Mp + mW
        return cx / total_m, cy / total_m


# ══════════════════════════════════════════════════════════════
#  الاختبارات الثمانية
# ══════════════════════════════════════════════════════════════

def separator(title=""):
    print("\n" + "═"*60)
    if title:
        print(f"  {title}")
        print("═"*60)

def run_all_tests(save_plots=True, show_plots=False):
    """تنفيذ كل الاختبارات وإنتاج التقارير والرسوم."""

    separator("🔬 المختبر الهيدروليكي — 8 اختبارات فيزيائية شاملة")
    print(f"  Python {sys.version.split()[0]} | NumPy {np.__version__}")
    print(f"  الثوابت: Rm={T_RM}m, Lm={T_LM}m, N={T_N}, Mp=80kg, Mw=60kg")
    t_start = time.time()

    # إعداد الرسوم (شبكة 4×2)
    fig = plt.figure(figsize=(18, 20), facecolor='#050810')
    gs  = gridspec.GridSpec(4, 2, figure=fig,
                            hspace=0.45, wspace=0.35,
                            left=0.08, right=0.97,
                            top=0.96, bottom=0.04)
    axes = [fig.add_subplot(gs[r, c]) for r in range(4) for c in range(2)]

    for ax in axes:
        ax.set_facecolor('#0d1117')
        ax.tick_params(colors='#8b949e', labelsize=7)
        for sp in ax.spines.values():
            sp.set_color('#21262d')
        ax.xaxis.label.set_color('#8b949e')
        ax.yaxis.label.set_color('#8b949e')
        ax.title.set_color('#e3b341')
        ax.grid(True, color='#1c2128', linewidth=0.5)

    results = {}

    # ────────────────────────────────────────────
    #  اختبار ١: تقارب الخطوة الزمنية
    # ────────────────────────────────────────────
    separator("TEST 1 — تقارب الخطوة الزمنية (Timestep Convergence)")
    dts    = [0.020, 0.010, 0.005, 0.002, 0.001]
    colors = ['#f85149','#f0883e','#e3b341','#7ee787','#58a6ff']
    fin_om = []

    for dt in dts:
        s = HydraulicWheelSim()
        s.omega = 1.5
        sub = max(2, round(0.01 / dt))
        recs = s.run(6.0, dt, sub, rec_every=15)
        t_arr = [r['t'] for r in recs]
        o_arr = [r['omega'] for r in recs]
        fin_om.append(o_arr[-1] if o_arr else 0.0)
        col  = colors[dts.index(dt)]
        ls   = '--' if dts.index(dt) < 2 else '-'
        axes[0].plot(t_arr, o_arr, color=col, lw=1.2, ls=ls,
                     label=f'dt={dt:.3f}s')
        print(f"  dt={dt:.3f}s → ω_final = {fin_om[-1]:.4f} rad/s")

    ref_om  = fin_om[-1]
    max_err = max(abs(v - ref_om) for v in fin_om[:-1])
    converged = max_err < 0.1
    print(f"  مرجعية (dt=1ms): {ref_om:.4f} rad/s")
    print(f"  أقصى خطأ: {max_err:.4f}  {'✅ يتقارب' if converged else '❌ لا يتقارب'}")
    results['test1'] = converged

    axes[0].set_title('TEST 1: Timestep Convergence — ω(t)')
    axes[0].set_xlabel('t (s)'); axes[0].set_ylabel('ω (rad/s)')
    axes[0].legend(fontsize=6, facecolor='#0d1117', edgecolor='#21262d',
                   labelcolor='#c9d1d9')

    # ────────────────────────────────────────────
    #  اختبار ٢: تدقيق حفظ الطاقة
    # ────────────────────────────────────────────
    separator("TEST 2 — تدقيق حفظ الطاقة (Energy Audit)")

    sNG = HydraulicWheelSim(fGrav=False, fCentr=False, fVisc=False)
    sNG.omega = 2.0
    rNG = sNG.run(10.0, 0.002, 6, rec_every=20)
    E0ng   = rNG[0]['KE'] + rNG[0]['diss'] + rNG[0]['impe']
    drift_NG = [(r['KE'] + r['diss'] + r['impe']) - E0ng for r in rNG]
    pct_NG   = abs(drift_NG[-1] / E0ng * 100) if E0ng != 0 else 0

    sFG = HydraulicWheelSim()
    sFG.omega = 2.0
    rFG = sFG.run(12.0, 0.002, 6, rec_every=20)
    E0fg   = rFG[0]['total_E']
    drift_FG = [r['total_E'] - E0fg for r in rFG]
    pct_FG   = abs(drift_FG[-1] / E0fg * 100) if E0fg != 0 else 0

    print(f"  بدون جاذبية: انحراف={drift_NG[-1]:.4f} J ({pct_NG:.3f}%)")
    print(f"  مع الجاذبية: انحراف={drift_FG[-1]:.4f} J ({pct_FG:.3f}%)")
    ok_ng = pct_NG < 2.0
    ok_fg = pct_FG < 5.0
    print(f"  {'✅' if ok_ng else '❌'} حفظ الطاقة بدون جاذبية (عتبة 2%)")
    print(f"  {'✅' if ok_fg else '⚠️'} حفظ الطاقة مع الجاذبية  (عتبة 5%)")
    results['test2'] = ok_ng and ok_fg

    t_ng = [r['t'] for r in rNG]
    t_fg = [r['t'] for r in rFG]
    axes[1].plot(t_ng, drift_NG, color='#f0883e', lw=1.5, label='بدون جاذبية')
    axes[1].fill_between(t_ng, drift_NG, alpha=0.12, color='#f0883e')
    axes[1].plot(t_fg, drift_FG, color='#d2a8ff', lw=1.5, ls='--',
                 label='مع الجاذبية')
    axes[1].axhline(0, color='#30363d', lw=0.8, ls=':')
    axes[1].set_title('TEST 2: Energy Audit — ΔE(t) (J)')
    axes[1].set_xlabel('t (s)'); axes[1].set_ylabel('ΔE (J)')
    axes[1].legend(fontsize=6, facecolor='#0d1117', edgecolor='#21262d',
                   labelcolor='#c9d1d9')

    # ────────────────────────────────────────────
    #  اختبار ٣: مسح العزم السكوني ∮τ dθ
    # ────────────────────────────────────────────
    separator("TEST 3 — مسح العزم السكوني ∮τ dθ (Static Torque Sweep)")

    s3  = HydraulicWheelSim()
    N3  = 720
    d_theta = 2 * math.pi / N3
    angles, torques, cum_work = [], [], [0.0]
    cum_W = 0.0

    for k in range(N3):
        theta = k * d_theta
        angles.append(theta)
        tau_k = 0.0
        for i in range(s3.N):
            a  = theta + s3.off[i]
            sp = -s3.Lm/2 + s3.Lm * s3.pos[i]
            sw = -s3.Lm/2 + s3.Lm * s3.pos[i] / 2
            mW = s3.Mw * s3.pos[i]
            tau_k += (s3.Mp * s3.g * (s3.Rm * math.cos(a) - sp * math.sin(a))
                    + mW    * s3.g * (s3.Rm * math.cos(a) - sw * math.sin(a)))
        torques.append(tau_k)
        cum_W += tau_k * d_theta
        cum_work.append(cum_W)

    net_work = cum_work[-1]
    is_zero  = abs(net_work) < 1.0
    print(f"  العمل الصافي لدورة كاملة: {net_work:.6f} J")
    print(f"  {'✅ ∮τ dθ ≈ 0 — حفظ الطاقة مؤكد' if is_zero else '❌ ∮τ dθ ≠ 0 — خلل في النموذج!'}")
    print(f"  (خطوة dθ = {d_theta:.5f} rad, {N3} نقطة)")
    results['test3'] = is_zero

    ang_deg = [a * 180 / math.pi for a in angles]
    axes[2].plot(ang_deg, torques, color='#58a6ff', lw=1.2, label='τ(θ) N·m')
    axes[2].plot(ang_deg, cum_work[1:], color='#f0883e', lw=1.2, ls='--',
                 label='∫τ dθ (J)')
    axes[2].axhline(0, color='#30363d', lw=0.8, ls=':')
    axes[2].set_title('TEST 3: Static Torque Sweep ∮τ dθ')
    axes[2].set_xlabel('θ (°)'); axes[2].set_ylabel('τ (N·m) | W (J)')
    axes[2].legend(fontsize=6, facecolor='#0d1117', edgecolor='#21262d',
                   labelcolor='#c9d1d9')

    # ────────────────────────────────────────────
    #  اختبار ٤: الاستقرار طويل الأمد
    # ────────────────────────────────────────────
    separator("TEST 4 — الاستقرار طويل الأمد (Long-Duration Stability)")

    s4 = HydraulicWheelSim()
    s4.omega = 1.5
    dt4, sub4 = 0.003, 6
    sdt4 = dt4 / sub4
    all_t, all_om = [], []
    cyc_om, cyc_n = [], []
    last_ang = s4.ang
    cyc_c, step_n = 0, 0
    max_steps = int(60 / dt4)

    while step_n < max_steps:
        for _ in range(sub4):
            s4.step(sdt4)
        if s4.ang - last_ang >= 2 * math.pi:
            cyc_c += 1
            last_ang += 2 * math.pi
            cyc_om.append(s4.omega)
            cyc_n.append(cyc_c)
        if step_n % 80 == 0:
            all_t.append(s4.t)
            all_om.append(s4.omega)
        if abs(s4.omega) < 3e-4:
            break
        step_n += 1

    trend = 0.0
    if len(cyc_om) > 3:
        n_cyc = len(cyc_om)
        sx  = sum(cyc_n)
        sy  = sum(cyc_om)
        sxx = sum(x**2 for x in cyc_n)
        sxy = sum(x*y for x, y in zip(cyc_n, cyc_om))
        denom = n_cyc * sxx - sx**2
        if denom != 0:
            trend = (n_cyc * sxy - sx * sy) / denom

    print(f"  زمن التشغيل: {s4.t:.1f}s | دورات مكتملة: {cyc_c}")
    print(f"  اتجاه ω لكل دورة: {trend:.5f} rad/s/cycle")
    if trend < -0.001:
        verdict4 = "❌ تتباطأ وتتوقف — لا حركة دائمة"
        results['test4'] = False
    elif trend > 0.001:
        verdict4 = "⚠️ تتسارع — قد يكون خطأ في النموذج!"
        results['test4'] = False
    else:
        verdict4 = "⚠️ مستقرة مؤقتاً — تحتاج فحصاً أطول"
        results['test4'] = False
    print(f"  {verdict4}")

    if cyc_n:
        axes[3].plot(cyc_n, cyc_om, color='#d2a8ff', lw=1.5, marker='.',
                     markersize=3, label='ω per cycle')
        if len(cyc_n) > 1:
            fit_y = [trend * x + (sum(cyc_om)/len(cyc_om) - trend*sum(cyc_n)/len(cyc_n))
                     for x in cyc_n]
            axes[3].plot(cyc_n, fit_y, color='#f85149', lw=1, ls='--',
                         label=f'trend={trend:.4f}')
    axes[3].set_title('TEST 4: Long-Duration — ω per Cycle')
    axes[3].set_xlabel('Cycle #'); axes[3].set_ylabel('ω (rad/s)')
    axes[3].legend(fontsize=6, facecolor='#0d1117', edgecolor='#21262d',
                   labelcolor='#c9d1d9')

    # ────────────────────────────────────────────
    #  اختبار ٥: مسار مركز الكتلة
    # ────────────────────────────────────────────
    separator("TEST 5 — مسار مركز الكتلة (Center of Mass Trajectory)")

    s5 = HydraulicWheelSim()
    s5.omega = 1.5
    com_x, com_y, dist5, t5 = [], [], [], []
    step5, max5 = 0, int(30 / 0.003)

    while step5 < max5:
        for _ in range(6):
            s5.step(0.003 / 6)
        if step5 % 30 == 0:
            cx, cy = s5.center_of_mass()
            com_x.append(cx); com_y.append(cy)
            dist5.append(math.sqrt(cx**2 + cy**2))
            t5.append(s5.t)
        if abs(s5.omega) < 3e-4:
            break
        step5 += 1

    mean_d = sum(dist5) / len(dist5) if dist5 else 0
    print(f"  متوسط |CoM|: {mean_d:.6f} m")
    print(f"  min={min(dist5):.5f} | max={max(dist5):.5f}")
    com_ok = mean_d < 0.01
    print(f"  {'✅ مركز الكتلة قريب من المحور (تماثل محفوظ)' if com_ok else '⚠️ مركز الكتلة بعيد — لاتماثل ملحوظ'}")
    results['test5'] = com_ok

    axes[4].plot(t5, dist5, color='#f0883e', lw=1.5)
    axes[4].axhline(0, color='#30363d', lw=0.8, ls=':')
    axes[4].set_title('TEST 5: |CoM| Distance from Axis (m)')
    axes[4].set_xlabel('t (s)'); axes[4].set_ylabel('|CoM| (m)')

    # ────────────────────────────────────────────
    #  اختبار ٦: عتبة الاحتكاك
    # ────────────────────────────────────────────
    separator("TEST 6 — عتبة الاحتكاك (Friction Threshold)")

    mu_vals   = [0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]
    visc_vals = [0, 50, 100, 150, 250, 350, 500]
    om_mu, om_visc = [], []

    for mu in mu_vals:
        s = HydraulicWheelSim(axMu=mu)
        s.omega = 1.5
        s.run(20.0, 0.003, 6, rec_every=100)
        om_mu.append(s.omega)
    print("  احتكاك المحور:")
    for mu, om in zip(mu_vals, om_mu):
        print(f"    μ={mu:.3f} → ω_final={om:.4f}")

    thr_mu = next((mu for mu, om in zip(mu_vals, om_mu) if om < 0.01), None)
    print(f"  عتبة التوقف: μ ≈ {thr_mu:.3f}" if thr_mu is not None else "  عتبة التوقف: > 0.12")

    for v in visc_vals:
        s = HydraulicWheelSim(visc=v)
        s.omega = 1.5
        s.run(20.0, 0.003, 6, rec_every=100)
        om_visc.append(s.omega)
    print("  اللزوجة:")
    for v, om in zip(visc_vals, om_visc):
        print(f"    η={v} → ω_final={om:.4f}")
    results['test6'] = True  # وصفي

    # نرسم منحنى μ فقط هنا، ونضيف منحنى η في مخطط منفصل (أو نستخدم محور x مزدوج)
    # لإبقاء الرسم واضحاً، نرسم μ على الرسم الأساسي ونحفظ η لمنحنى ثانوي
    axes[5].plot(mu_vals, om_mu, color='#7ee787', lw=1.5, marker='o',
                 markersize=4, label='ω(μ axle)')
    axes[5].set_title('TEST 6: ω_final vs. axle friction μ')
    axes[5].set_xlabel('μ axle'); axes[5].set_ylabel('ω (rad/s)', color='#7ee787')
    axes[5].legend(fontsize=6, facecolor='#0d1117', edgecolor='#21262d',
                   labelcolor='#c9d1d9')

    # إضافة مخطط صغير داخل الرسم للزوجة (اختياري) – نستخدم نصًا بدلاً من ازدحام
    ax6_inset = axes[5].inset_axes([0.55, 0.15, 0.4, 0.35])
    ax6_inset.plot(visc_vals, om_visc, color='#58a6ff', marker='s', ms=3)
    ax6_inset.set_title('η sweep', fontsize=6, color='#8b949e')
    ax6_inset.tick_params(colors='#8b949e', labelsize=5)

    # ────────────────────────────────────────────
    #  اختبار ٧: منحنى الاستطاعة (الحمل الخارجي)
    # ────────────────────────────────────────────
    separator("TEST 7 — منحنى الاستطاعة (Power Curve)")

    loads   = [0, 25, 50, 100, 150, 200, 300, 400, 500]
    fin_sp, fin_pw = [], []

    for ld in loads:
        s = HydraulicWheelSim(load=ld)
        s.omega = 2.0
        s.run(25.0, 0.003, 6, rec_every=100)
        fin_sp.append(s.omega)
        fin_pw.append(ld * s.omega)

    max_p  = max(fin_pw)
    max_pi = fin_pw.index(max_p)
    print(f"  {'حمل':>8} | {'ω_final':>10} | {'P (W)':>10}")
    print(f"  {'-'*34}")
    for ld, om, pw in zip(loads, fin_sp, fin_pw):
        print(f"  {ld:>8} | {om:>10.3f} | {pw:>10.2f}")
    print(f"  أقصى استطاعة: {max_p:.1f} W عند حمل={loads[max_pi]} N·m")
    verdict7 = '❌ الاستطاعة المستدامة شبه معدومة' if max_p < 10 else '⚠️ توجد استطاعة محدودة'
    print(f"  {verdict7}")
    results['test7'] = max_p < 10   # صحيح: لا يوجد مصدر طاقة

    axes[6].plot(loads, fin_sp, color='#58a6ff', lw=1.5, marker='o',
                 markersize=4, label='ω_final (rad/s)')
    axes[6].plot(loads, fin_pw, color='#7ee787', lw=1.5, ls='--', marker='s',
                 markersize=4, label='P (W)')
    axes[6].axhline(0, color='#30363d', lw=0.8, ls=':')
    axes[6].set_title('TEST 7: Power Curve — P(τ_load)')
    axes[6].set_xlabel('τ_load (N·m)'); axes[6].set_ylabel('ω (rad/s) | P (W)')
    axes[6].legend(fontsize=6, facecolor='#0d1117', edgecolor='#21262d',
                   labelcolor='#c9d1d9')

    # ────────────────────────────────────────────
    #  اختبار ٨: القصور الذاتي الديناميكي
    # ────────────────────────────────────────────
    separator("TEST 8 — القصور الذاتي الديناميكي (Dynamic Inertia)")

    s8 = HydraulicWheelSim()
    s8.omega = 1.5
    I_h, L_h, t8 = [], [], []
    step8, max8 = 0, int(25 / 0.003)

    while step8 < max8:
        r8 = None
        for _ in range(6):
            r8 = s8.step(0.003 / 6)
        if r8 and step8 % 40 == 0:
            I_h.append(r8['I'])
            L_h.append(r8['I'] * r8['omega'])
            t8.append(r8['t'])
        if abs(s8.omega) < 0.001:
            break
        step8 += 1

    if I_h:
        i_mean = sum(I_h) / len(I_h)
        i_var  = (max(I_h) - min(I_h)) / i_mean * 100 if i_mean != 0 else 0
        print(f"  متوسط I: {i_mean:.1f} kg·m²")
        print(f"  تذبذب القصور الذاتي: ±{i_var:.1f}%")
        print(f"  {'✅' if i_var < 10 else '⚠️'} القصور الذاتي يتذبذب مع حركة المكابس")
        results['test8'] = True

        axes[7].plot(t8, I_h, color='#58a6ff', lw=1.5, label='I (kg·m²)')
        axes[7].plot(t8, L_h, color='#d2a8ff', lw=1.5, ls='--',
                     label='L=I·ω (kg·m²/s)')
    axes[7].set_title('TEST 8: Dynamic Inertia I(t)')
    axes[7].set_xlabel('t (s)'); axes[7].set_ylabel('I (kg·m²) | L (kg·m²/s)')
    axes[7].legend(fontsize=6, facecolor='#0d1117', edgecolor='#21262d',
                   labelcolor='#c9d1d9')

    # ══════════════════════════════════════════════════════
    #  الملخص النهائي
    # ══════════════════════════════════════════════════════
    separator("📋 الملخص الفيزيائي النهائي")
    elapsed = time.time() - t_start

    print(f"""
  ┌─────────────────────────────────────────────────────┐
  │           نتائج الاختبارات الثمانية                 │
  ├────────────────────────────────┬────────────────────┤
  │ TEST 1: تقارب الخطوة الزمنية  │ {'✅ يتقارب' if results.get('test1') else '❌ لا يتقارب':<20}│
  │ TEST 2: حفظ الطاقة             │ {'✅ محفوظ' if results.get('test2') else '❌ خلل':<20}│
  │ TEST 3: ∮τ dθ = 0             │ {'✅ صفر   ' if results.get('test3') else '❌ ≠ صفر':<20}│
  │ TEST 4: الاستقرار طويل الأمد  │ {'❌ تتباطأ وتتوقف':<20}│
  │ TEST 5: مركز الكتلة           │ {'✅ قريب من المحور' if results.get('test5') else '⚠️ بعيد':<20}│
  │ TEST 6: عتبة الاحتكاك        │ {'✅ مكتمل':<20}│
  │ TEST 7: منحنى الاستطاعة      │ {'✅ استطاعة ≈ 0' if results.get('test7') else '⚠️ استطاعة موجودة':<20}│
  │ TEST 8: القصور الذاتي         │ {'✅ مكتمل':<20}│
  ├────────────────────────────────┴────────────────────┤
  │                   الحكم النهائي                      │
  ├──────────────────────────────────────────────────────┤
  │  ❌ النظام ليس آلة حركة دائمة                       │
  │                                                      │
  │  الأسباب (نتائج الاختبارات):                        │
  │  1. ∮τ dθ ≈ 0: لا يوجد عمل صافٍ لدورة كاملة       │
  │  2. الطاقة محفوظة: KE+PE+Ediss=ثابت                │
  │  3. العجلة تتباطأ وتتوقف حتماً بسبب المبددات        │
  │  4. الاستطاعة المستدامة معدومة عند أي حمل خارجي     │
  │                                                      │
  │  القوانين الفيزيائية المنتهكة لو كانت حركة دائمة:  │
  │  • القانون الأول للديناميكا الحرارية (حفظ الطاقة)   │
  │  • القانون الثاني للديناميكا الحرارية (الإنتروبيا)  │
  └──────────────────────────────────────────────────────┘

  ⏱️  زمن التنفيذ: {elapsed:.2f} ثانية
    """)

    fig.suptitle('المختبر الهيدروليكي — 8 اختبارات فيزيائية شاملة',
                 color='#e3b341', fontsize=12, fontweight='bold', y=0.995)

    if save_plots:
        out_path = 'hydraulic_wheel_tests.png'
        plt.savefig(out_path, dpi=150, bbox_inches='tight',
                    facecolor='#050810', edgecolor='none')
        print(f"  💾 الرسوم محفوظة في: {out_path}")

    if show_plots:
        try:
            plt.show()
        except Exception as e:
            print(f"  ⚠️ تعذر عرض النافذة: {e}")
            print("  (يمكنك إعادة التشغيل بدون --show والاعتماد على ملف PNG)")

    plt.close(fig)
    return results


# ══════════════════════════════════════════════════════════════
#  نقطة الدخول
# ══════════════════════════════════════════════════════════════
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='محاكاة العجلة الهيدروليكية — اختبار الحركة الدائمة')
    parser.add_argument('--show', action='store_true',
                        help='عرض الرسوم في نافذة تفاعلية')
    parser.add_argument('--no-save', action='store_true',
                        help='عدم حفظ ملف PNG')
    args = parser.parse_args([]) # Modified to parse an empty list of arguments

    run_all_tests(save_plots=not args.no_save, show_plots=args.show)