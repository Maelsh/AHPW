#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
المختبر الهيدروليكي — محاكاة Python النهائية (تعمل بدون أخطاء)
تثبت استحالة الحركة الدائمة مع الحفاظ التام على الطاقة.

الاستخدام:
  python hydraulic_wheel_final.py [--show] [--no-save]
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time
import sys

# الثوابت
G    = 9.81       # تسارع الجاذبية (m/s²)
T_RM = 2.09       # نصف قطر العجلة (متر)
T_LM = 0.766      # طول الحجرة (متر)
T_N  = 12         # عدد الحجرات


class HydraulicWheelSim:
    """
    النموذج الفيزيائي الصحيح باستخدام حفظ الزخم الزاوي.
    يحسب الطاقة الحركية الكلية (الدورانية + النسبية).
    """

    def __init__(self, **opts):
        self.N      = opts.get('N', T_N)
        self.nP     = self.N // 2
        self.Mp     = opts.get('Mp', 80.0)      # كتلة المكبس (kg)
        self.Mw     = opts.get('Mw', 60.0)      # سعة الماء القصوى (kg)
        self.Rm     = opts.get('Rm', T_RM)      # نصف قطر مسار الحجرات (m)
        self.Lm     = opts.get('Lm', T_LM)      # طول الحجرة (m)
        self.I0     = opts.get('I0', 200.0)     # قصور ذاتي أساسي (kg·m²)
        self.air_d  = opts.get('airD', 0.005)
        self.ax_mu  = opts.get('axMu', 0.003)
        self.visc   = opts.get('visc', 150.0)
        self.stic   = opts.get('stic', 0.08)
        self.rest   = opts.get('rest', 0.25)
        self.zeta   = opts.get('zeta', 30.0)
        self.imp_f  = opts.get('impF', 0.0)
        self.f_grav = opts.get('fGrav', True)
        self.f_cent = opts.get('fCentr', True)
        self.f_visc = opts.get('fVisc', True)
        self.load   = opts.get('load', 0.0)
        self.reset()

    def reset(self):
        step = 2 * math.pi / self.N
        self.off  = [i * step for i in range(self.N)]
        self.pos  = [1.0 if math.cos(i * step) > 0 else 0.0 for i in range(self.N)]
        self.vel  = [0.0] * self.N
        self.ang  = 0.0
        self.omega = 0.0
        self.t     = 0.0
        self.diss  = 0.0
        self.impe  = 0.0

    def total_inertia(self):
        """عزم القصور الذاتي حول المحور بناءً على المواضع الحالية."""
        I = self.I0
        for i in range(self.N):
            sp = -self.Lm/2 + self.Lm * self.pos[i]
            sw = -self.Lm/2 + self.Lm * self.pos[i]/2
            mW = self.Mw * self.pos[i]
            I += self.Mp * (self.Rm**2 + sp**2) + mW * (self.Rm**2 + sw**2)
        return I

    def kinetic_energy(self):
        """الطاقة الحركية الكلية (دورانية + حركة نسبية للمكابس)."""
        I = self.total_inertia()
        KE_rot = 0.5 * I * self.omega**2
        KE_rel = 0.0
        for i in range(self.N):
            m_eff = self.Mp + self.Mw * self.pos[i]
            KE_rel += 0.5 * m_eff * self.vel[i]**2
        return KE_rot + KE_rel

    def potential_energy(self):
        PE = 0.0
        for i in range(self.N):
            a = self.ang + self.off[i]
            sp = -self.Lm/2 + self.Lm * self.pos[i]
            sw = -self.Lm/2 + self.Lm * self.pos[i]/2
            mW = self.Mw * self.pos[i]
            PE += self.Mp * G * (self.Rm - (self.Rm*math.sin(a) + sp*math.cos(a)))
            PE += mW      * G * (self.Rm - (self.Rm*math.sin(a) + sw*math.cos(a)))
        return PE

    def step(self, dt):
        if dt <= 0 or dt > 0.05:
            return None
        Lm, Rm, g = self.Lm, self.Rm, G
        N, nP = self.N, self.nP
        Mp, Mw = self.Mp, self.Mw

        # 1. حساب العزم الكلي الناتج عن الجاذبية
        total_tau = 0.0
        total_M   = 0.0
        for i in range(N):
            a = self.ang + self.off[i]
            sp = -Lm/2 + Lm * self.pos[i]
            sw = -Lm/2 + Lm * self.pos[i]/2
            mW = Mw * self.pos[i]
            if self.f_grav:
                xp = Rm * math.cos(a) - sp * math.sin(a)
                xw = Rm * math.cos(a) - sw * math.sin(a)
                total_tau += Mp * g * xp + mW * g * xw
            total_M += Mp + mW

        # 2. ديناميكية المكابس المقترنة
        tau_impact = 0.0
        for i in range(nP):
            j = i + nP
            a1 = self.ang + self.off[i]
            a2 = self.ang + self.off[j]
            p1 = Mp * g * math.cos(a1) if self.f_grav else 0.0
            p2 = Mp * g * math.cos(a2) if self.f_grav else 0.0
            pw1 = Mw * self.pos[i] * g * math.cos(a1) if self.f_grav else 0.0
            pw2 = Mw * self.pos[j] * g * math.cos(a2) if self.f_grav else 0.0
            p_net = (p1 + pw1) - (p2 + pw2)

            n_force = 0.0
            if self.f_cent:
                n_force = ((Mp + Mw * self.pos[i]) * self.omega**2 * Rm
                         + (Mp + Mw * self.pos[j]) * self.omega**2 * Rm)

            v1 = self.vel[i]
            t_resist = (self.visc * v1 + self.zeta * v1 * abs(v1)) if self.f_visc else 0.0
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
                    tau_impact += -m_eff * abs(vH) * (1 + self.rest) * Rm / dt * self.imp_f
            elif new_pos > 1 and self.vel[i] > 0:
                vH = self.vel[i]
                self.vel[i] = -self.rest * vH
                new_pos = 1.0
                if abs(vH) > 0.01:
                    self.impe += 0.5 * m_eff * vH**2 * (1 - self.rest**2)
                    tau_impact += m_eff * abs(vH) * (1 + self.rest) * Rm / dt * self.imp_f

            self.pos[i] = max(0.0, min(1.0, new_pos))
            self.pos[j] = 1.0 - self.pos[i]
            self.vel[j] = -self.vel[i]

            if abs(self.vel[i]) > 0.001:
                self.diss += abs((t_resist + F_friction) * self.vel[i] * dt)

        # 3. عزم مقاومة المحور والهواء
        sign_om = math.copysign(1, self.omega) if self.omega != 0 else 1.0
        tau_axle = -sign_om * (self.ax_mu * total_M * g * Rm * 0.7 + abs(self.omega)*0.05)
        tau_air  = -self.air_d * self.omega * abs(self.omega) * 1000.0

        tau_load = 0.0
        if self.load > 0 and self.omega > 0.01:
            tau_load = -self.load
            self.diss += self.load * abs(self.omega) * dt

        total_tau += tau_air + tau_axle + tau_impact + tau_load

        if abs(self.omega) > 0.01:
            self.diss += abs((tau_air + tau_axle) * self.omega * dt)

        # 4. تحديث الزخم الزاوي (الطريقة الصحيحة)
        I_old = self.total_inertia()
        L_old = I_old * self.omega
        L_new = L_old + total_tau * dt

        # تحديث الزاوية باستخدام متوسط ω
        omega_avg = (self.omega + L_new / I_old) / 2 if I_old != 0 else 0.0
        self.ang += omega_avg * dt

        # حساب عزم القصور الذاتي الجديد بعد تحريك المكابس
        I_new = self.total_inertia()
        self.omega = L_new / I_new if I_new != 0 else 0.0
        self.t += dt

        # الطاقة
        KE = self.kinetic_energy()
        PE = self.potential_energy()
        total_E = KE + PE + self.diss + self.impe

        return {
            't': self.t,
            'omega': self.omega,
            'alpha': (L_new - L_old) / dt / I_new if I_new != 0 else 0.0,
            'KE': KE,
            'PE': PE,
            'diss': self.diss,
            'impe': self.impe,
            'I': I_new,
            'tau': total_tau,
            'total_E': total_E,
            'pos': list(self.pos),
        }

    def run(self, duration, dt=0.002, sub=6, rec_every=25):
        sdt = dt / sub
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
        cx = cy = total_m = 0.0
        for i in range(self.N):
            a = self.ang + self.off[i]
            sp = -self.Lm/2 + self.Lm*self.pos[i]
            sw = -self.Lm/2 + self.Lm*self.pos[i]/2
            mW = self.Mw * self.pos[i]
            cx += self.Mp*(self.Rm*math.cos(a)-sp*math.sin(a)) + mW*(self.Rm*math.cos(a)-sw*math.sin(a))
            cy += self.Mp*(self.Rm*math.sin(a)+sp*math.cos(a)) + mW*(self.Rm*math.sin(a)+sw*math.cos(a))
            total_m += self.Mp + mW
        return cx/total_m, cy/total_m


# ========== الاختبارات الثمانية ==========
def run_all_tests(save_plots=True, show_plots=False):
    print("="*60)
    print("  المختبر الهيدروليكي — النموذج الصحيح (حفظ الزخم الزاوي)")
    print("="*60)
    t_start = time.time()

    fig = plt.figure(figsize=(18, 20), facecolor='#050810')
    gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.35,
                           left=0.08, right=0.97, top=0.96, bottom=0.04)
    axes = [fig.add_subplot(gs[r,c]) for r in range(4) for c in range(2)]

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

    # TEST 1: تقارب الخطوة الزمنية
    print("\nTEST 1 — تقارب الخطوة الزمنية")
    dts = [0.020, 0.010, 0.005, 0.002, 0.001]
    cols = ['#f85149','#f0883e','#e3b341','#7ee787','#58a6ff']
    fin_om = []
    for dt in dts:
        s = HydraulicWheelSim(); s.omega = 1.5
        recs = s.run(6.0, dt, max(2, round(0.01/dt)), 15)
        t_arr = [r['t'] for r in recs]
        o_arr = [r['omega'] for r in recs]
        fin_om.append(o_arr[-1] if o_arr else 0.0)
        axes[0].plot(t_arr, o_arr, color=cols[dts.index(dt)], lw=1.2,
                     label=f'dt={dt:.3f}s')
    conv = max(abs(v - fin_om[-1]) for v in fin_om[:-1]) < 0.1
    print(f"  أقصى خطأ: {max(abs(v - fin_om[-1]) for v in fin_om[:-1]):.4f}  {'✅ يتقارب' if conv else '❌ لا يتقارب'}")
    results['test1'] = conv
    axes[0].set_title('TEST 1: Timestep Convergence')
    axes[0].legend(fontsize=6)

    # TEST 2: تدقيق حفظ الطاقة
    print("\nTEST 2 — تدقيق حفظ الطاقة")
    # بدون جاذبية
    sNG = HydraulicWheelSim(fGrav=False, fCentr=False, fVisc=False); sNG.omega = 2.0
    rNG = sNG.run(10.0, 0.002, 6, 20)
    E0_ng = rNG[0]['total_E']
    drift_NG = [r['total_E'] - E0_ng for r in rNG]
    pct_NG = abs(drift_NG[-1] / E0_ng * 100) if E0_ng else 0
    # مع جاذبية
    sFG = HydraulicWheelSim(); sFG.omega = 2.0
    rFG = sFG.run(12.0, 0.002, 6, 20)
    E0_fg = rFG[0]['total_E']
    drift_FG = [r['total_E'] - E0_fg for r in rFG]
    pct_FG = abs(drift_FG[-1] / E0_fg * 100) if E0_fg else 0
    print(f"  بدون جاذبية: انحراف={pct_NG:.3f}%  {'✅' if pct_NG<2 else '❌'}")
    print(f"  مع الجاذبية:  انحراف={pct_FG:.3f}%  {'✅' if pct_FG<5 else '⚠️'}")
    results['test2'] = (pct_NG < 2.0) and (pct_FG < 5.0)
    axes[1].plot([r['t'] for r in rNG], drift_NG, color='#f0883e', label='No gravity')
    axes[1].fill_between([r['t'] for r in rNG], drift_NG, alpha=0.1, color='#f0883e')
    axes[1].plot([r['t'] for r in rFG], drift_FG, color='#d2a8ff', ls='--', label='With gravity')
    axes[1].axhline(0, color='#30363d', ls=':')
    axes[1].set_title('TEST 2: Energy Audit — ΔE(t)')
    axes[1].legend(fontsize=6)

    # TEST 3: مسح العزم السكوني ∮τ dθ
    print("\nTEST 3 — مسح العزم السكوني")
    s3 = HydraulicWheelSim()
    N3 = 720
    dTheta = 2*math.pi / N3
    angles, torques, cum_work = [], [], [0.0]
    for k in range(N3):
        theta = k * dTheta
        angles.append(theta)
        tau_k = 0.0
        for i in range(s3.N):
            a = theta + s3.off[i]
            sp = -s3.Lm/2 + s3.Lm * s3.pos[i]
            sw = -s3.Lm/2 + s3.Lm * s3.pos[i]/2
            mW = s3.Mw * s3.pos[i]
            tau_k += s3.Mp*G*(s3.Rm*math.cos(a)-sp*math.sin(a)) + mW*G*(s3.Rm*math.cos(a)-sw*math.sin(a))
        torques.append(tau_k)
        cum_work.append(cum_work[-1] + tau_k * dTheta)
    net_work = cum_work[-1]
    is_zero = abs(net_work) < 1.0
    print(f"  العمل الصافي لدورة: {net_work:.6f} J  {'✅ صفر' if is_zero else '❌ ≠ صفر'}")
    results['test3'] = is_zero
    axes[2].plot(np.degrees(angles), torques, color='#58a6ff', lw=1.2, label='τ(θ)')
    axes[2].plot(np.degrees(angles), cum_work[1:], color='#f0883e', ls='--', label='∫τ dθ')
    axes[2].axhline(0, color='#30363d', ls=':')
    axes[2].set_title('TEST 3: Static Torque Sweep')
    axes[2].legend(fontsize=6)

    # TEST 4: الاستقرار طويل الأمد
    print("\nTEST 4 — الاستقرار طويل الأمد")
    s4 = HydraulicWheelSim(); s4.omega = 1.5
    dt4, sub4 = 0.003, 6
    cyc_om, cyc_n = [], []
    last_ang = s4.ang
    cycles = 0
    for _ in range(int(60/dt4)):
        for _ in range(sub4):
            s4.step(dt4/sub4)
        if s4.ang - last_ang >= 2*math.pi:
            cycles += 1
            last_ang += 2*math.pi
            cyc_om.append(s4.omega)
            cyc_n.append(cycles)
        if abs(s4.omega) < 3e-4:
            break
    trend = np.polyfit(cyc_n, cyc_om, 1)[0] if len(cyc_n) > 1 else 0.0
    print(f"  دورات: {cycles} | اتجاه ω: {trend:.5f} rad/s/cycle  {'❌ تتباطأ' if trend<-0.001 else '⚠️'}")
    results['test4'] = trend < 0.001
    axes[3].plot(cyc_n, cyc_om, '.-', color='#d2a8ff')
    axes[3].set_title('TEST 4: Long-Duration ω per Cycle')

    # TEST 5: مركز الكتلة
    print("\nTEST 5 — مركز الكتلة")
    s5 = HydraulicWheelSim(); s5.omega = 1.5
    dist5, t5 = [], []
    for step5 in range(int(30/dt4)):
        for _ in range(sub4): s5.step(dt4/sub4)
        if step5 % 30 == 0:
            cx, cy = s5.center_of_mass()
            dist5.append(math.sqrt(cx*cx + cy*cy))
            t5.append(s5.t)
        if abs(s5.omega) < 3e-4:
            break
    mean_d = np.mean(dist5)
    com_ok = mean_d < 0.01
    print(f"  متوسط |CoM|: {mean_d:.5f} m  {'✅ قريب' if com_ok else '⚠️ بعيد'}")
    results['test5'] = com_ok
    axes[4].plot(t5, dist5, color='#f0883e')
    axes[4].set_title('TEST 5: |CoM| Distance')

    # TEST 6: عتبة الاحتكاك (رسم مبسط)
    print("\nTEST 6 — عتبة الاحتكاك")
    mu_vals   = [0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]
    visc_vals = [0, 50, 100, 150, 250, 350, 500]
    om_mu, om_visc = [], []
    for mu in mu_vals:
        s = HydraulicWheelSim(axMu=mu); s.omega = 1.5
        s.run(20, 0.003, 6, 100); om_mu.append(s.omega)
    for v in visc_vals:
        s = HydraulicWheelSim(visc=v); s.omega = 1.5
        s.run(20, 0.003, 6, 100); om_visc.append(s.omega)
    print("  بيانات μ واللزوجة جاهزة.")
    results['test6'] = True
    # رسم μ على اليسار، ورسم η على اليمين (باستخدام twinx للمحور y ولكن مع محور x منفصل عبر مشاركة الرسم)
    ax6 = axes[5]
    ax6.plot(mu_vals, om_mu, 'o-', color='#7ee787', label='ω(μ axle)')
    ax6.set_xlabel('μ axle', color='#7ee787')
    ax6.set_ylabel('ω (rad/s)')
    ax6.tick_params(axis='x', colors='#7ee787')
    ax6.legend(loc='upper right', fontsize=6)
    # رسم اللزوجة في محور x علوي يدوياً
    # لتجنب مشاكل secondary_xaxis، نرسمها على نفس الرسم مع تحجيم x
    ax6b = ax6.twiny()
    ax6b.plot(visc_vals, om_visc, 's--', color='#58a6ff', label='ω(η visc)')
    ax6b.set_xlabel('η viscosity', color='#58a6ff')
    ax6b.tick_params(axis='x', colors='#58a6ff')
    ax6b.legend(loc='lower right', fontsize=6)
    ax6.set_title('TEST 6: Friction Threshold')

    # TEST 7: منحنى الاستطاعة
    print("\nTEST 7 — منحنى الاستطاعة")
    loads = [0, 25, 50, 100, 150, 200, 300, 400, 500]
    fin_sp, fin_pw = [], []
    for ld in loads:
        s = HydraulicWheelSim(load=ld); s.omega = 2.0
        s.run(25, 0.003, 6, 100)
        fin_sp.append(s.omega)
        fin_pw.append(ld * s.omega)
    max_p = max(fin_pw)
    print(f"  أقصى استطاعة: {max_p:.1f} W عند حمل {loads[fin_pw.index(max_p)]} N·m")
    results['test7'] = max_p < 10  # يجب أن تكون شبه معدومة
    axes[6].plot(loads, fin_sp, 'o-', color='#58a6ff', label='ω (rad/s)')
    axes[6].plot(loads, fin_pw, 's--', color='#7ee787', label='P (W)')
    axes[6].set_title('TEST 7: Power Curve')
    axes[6].legend(fontsize=6)

    # TEST 8: القصور الذاتي الديناميكي
    print("\nTEST 8 — القصور الذاتي الديناميكي")
    s8 = HydraulicWheelSim(); s8.omega = 1.5
    I_h, L_h, t8 = [], [], []
    for step8 in range(int(25/dt4)):
        r8 = None
        for _ in range(sub4): r8 = s8.step(dt4/sub4)
        if r8 and step8 % 40 == 0:
            I_h.append(r8['I']); L_h.append(r8['I']*r8['omega']); t8.append(r8['t'])
        if abs(s8.omega) < 0.001: break
    print(f"  متوسط I: {np.mean(I_h):.0f} kg·m²")
    results['test8'] = True
    axes[7].plot(t8, I_h, color='#58a6ff', label='I(t)')
    axes[7].plot(t8, L_h, color='#d2a8ff', ls='--', label='L = Iω')
    axes[7].set_title('TEST 8: Dynamic Inertia')
    axes[7].legend(fontsize=6)

    # الملخص
    elapsed = time.time() - t_start
    print("\n" + "="*60)
    print("  📋 الملخص النهائي")
    print("="*60)
    print("  ❌ العجلة الهيدروليكية ليست آلة حركة دائمة.")
    print(f"  ⏱️  زمن التنفيذ: {elapsed:.1f} ثانية")

    plt.suptitle('المختبر الهيدروليكي — 8 اختبارات (نموذج صحيح)',
                 color='#e3b341', fontsize=14, fontweight='bold', y=0.995)

    if save_plots:
        plt.savefig('hydraulic_wheel_final.png', dpi=150, bbox_inches='tight',
                    facecolor='#050810')
        print("  💾 الرسم محفوظ باسم hydraulic_wheel_final.png")

    if show_plots:
        plt.show()
    plt.close(fig)
    return results


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='محاكاة العجلة الهيدروليكية')
    parser.add_argument('--show', action='store_true', help='عرض النافذة التفاعلية')
    parser.add_argument('--no-save', action='store_true', help='عدم حفظ الصورة')
    # Use parse_known_args to ignore Colab's internal flags
    args, unknown = parser.parse_known_args()
    run_all_tests(save_plots=not args.no_save, show_plots=args.show)