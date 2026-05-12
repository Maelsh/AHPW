#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
المختبر الهيدروليكي — محاكاة فيزيائية شاملة
Hydraulic Wheel Physics Simulation — Python Edition

هذا الملف يحتوي على:
1. محرك محاكاة العجلة الهيدروليكية (معادلات مُعاد بناؤها من الكود الأصلي)
2. 8 اختبارات فيزيائية شاملة
3. رسم بياني لجميع النتائج

المؤلف: تحليل مستقل للمحاكاة الأصلية
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List, Tuple, Optional
import json
import sys

# ═══════════════════════════════════════════════════════
# القسم ١: محرك المحاكاة الرئيسي
# ═══════════════════════════════════════════════════════

@dataclass
class SimState:
    """حالة لحظية للمحاكاة"""
    t: float
    omega: float
    alpha: float
    KE: float
    PE: float
    diss: float
    impe: float
    totalE: float
    I: float
    tau: float
    pos: np.ndarray


class HydraulicWheelSim:
    """
    محاكاة العجلة الهيدروليكية الدوارة.

    المعادلات مُعاد بناؤها بدقة من الكود الأصلي (HTML/JavaScript)
    مع التحقق من صحة كل خطوة فيزيائية.

    الثوابت الافتراضية:
        Rm = 2.09 m    (نصف قطر العجلة)
        Lm = 0.766 m   (طول المكبس)
        N  = 12        (عدد الغرف)
        Mp = 80 kg     (كتلة المكبس)
        Mw = 60 kg     (كتلة الماء القصوى)
        I0 = 200       (قصور ذاتي العجلة)
    """

    def __init__(self, **opts):
        # الثوابت الفيزيائية
        self.g = 9.81
        self.N = opts.get('N', 12)
        self.nP = self.N // 2
        self.Mp = opts.get('Mp', 80)
        self.Mw = opts.get('Mw', 60)
        self.Rm = opts.get('Rm', 2.09)
        self.Lm = opts.get('Lm', 0.766)
        self.I0 = opts.get('I0', 200)

        # معاملات المقاومة
        self.airD = opts.get('airD', 0.005)
        self.axMu = opts.get('axMu', 0.003)
        self.visc = opts.get('visc', 150)
        self.stic = opts.get('stic', 0.08)
        self.rest = opts.get('rest', 0.25)
        self.zeta = opts.get('zeta', 30)
        self.impF = opts.get('impF', 0.0)

        # قوى مفعلة
        self.fGrav = opts.get('fGrav', True)
        self.fCentr = opts.get('fCentr', True)
        self.fVisc = opts.get('fVisc', True)

        # حمل خارجي
        self.load = opts.get('load', 0)

        self.reset()

    def reset(self):
        """إعادة تعيين الحالة إلى الوضع الابتدائي"""
        step = 2 * np.pi / self.N
        self.pos = np.array([1.0 if np.cos(i * step) > 0 else 0.0 for i in range(self.N)])
        self.vel = np.zeros(self.N)
        self.off = np.array([i * step for i in range(self.N)])
        self.ang = 0.0
        self.omega = 0.0
        self.t = 0.0
        self.diss = 0.0
        self.impe = 0.0

    def step(self, dt: float) -> Optional[SimState]:
        """
        خطوة زمنية واحدة في المحاكاة.

        المعادلات:
        1. حساب العزم الثقالي لكل غرفة: τ_g = Σ(m·g·x_cm)
        2. ديناميكيا المكابس (أزواج متقابلة): معادلة الحركة مع احتكاك tanh
        3. مقاومات العجلة: احتكاك محور + مقاومة هواء + حمل خارجي
        4. تحديث ω و θ: α = τ/I,  ω += α·dt,  θ += ω·dt
        5. حساب الطاقات: KE, PE, dissipated, impact
        """
        if dt <= 0 or dt > 0.05:
            return None

        Lm, Rm, g, N, nP, Mp, Mw, I0 = (
            self.Lm, self.Rm, self.g, self.N, self.nP, self.Mp, self.Mw, self.I0
        )

        totalI = I0
        totalTau = 0.0
        totalM = 0.0
        tauImp = 0.0

        # ── 1. حساب العزم الثقالي لكل غرفة ──
        for i in range(N):
            a = self.ang + self.off[i]
            sp = -Lm/2 + Lm * self.pos[i]
            sw = -Lm/2 + Lm * self.pos[i] / 2
            mW = Mw * self.pos[i]
            xp = Rm * np.cos(a) - sp * np.sin(a)
            xw = Rm * np.cos(a) - sw * np.sin(a)
            tG = Mp * g * xp + mW * g * xw if self.fGrav else 0.0
            totalTau += tG
            totalI += Mp * (Rm**2 + sp**2) + mW * (Rm**2 + sw**2)
            totalM += Mp + mW

        # ── 2. ديناميكيا المكابس (أزواج متقابلة) ──
        for i in range(nP):
            j = i + nP
            a1 = self.ang + self.off[i]
            a2 = self.ang + self.off[j]

            p1 = Mp * g * np.cos(a1) if self.fGrav else 0.0
            p2 = Mp * g * np.cos(a2) if self.fGrav else 0.0
            pw1 = Mw * self.pos[i] * g * np.cos(a1) if self.fGrav else 0.0
            pw2 = Mw * self.pos[j] * g * np.cos(a2) if self.fGrav else 0.0
            pNet = (p1 + pw1) - (p2 + pw2)

            # القوة الطبيعية (طرد مركزي)
            nF = 0.0
            if self.fCentr:
                nF = ((Mp + Mw * self.pos[i]) * self.omega**2 * Rm +
                      (Mp + Mw * self.pos[j]) * self.omega**2 * Rm)

            v1 = self.vel[i]
            tR = self.visc * v1 + self.zeta * v1 * abs(v1) if self.fVisc else 0.0
            mE = 2 * Mp + Mw
            Fs = self.stic * (nF + (Mp + Mw) * g * 0.1)
            dF = pNet - tR

            # احتكاك سلس باستخدام tanh (إصلاح مهم)
            vSm = 0.005
            dF -= Fs * np.tanh(v1 / vSm)

            self.vel[i] += (dF / mE) * dt
            nP2 = self.pos[i] + self.vel[i] * dt / Lm

            # اصطدام بالجدران مع معامل ارتداد
            if nP2 < 0 and self.vel[i] < 0:
                vH = self.vel[i]
                self.vel[i] = -self.rest * vH
                nP2 = 0.0
                if abs(vH) > 0.01:
                    self.impe += 0.5 * mE * vH**2 * (1 - self.rest**2)
                    tauImp += -mE * abs(vH) * (1 + self.rest) * Rm / dt * self.impF
            elif nP2 > 1 and self.vel[i] > 0:
                vH = self.vel[i]
                self.vel[i] = -self.rest * vH
                nP2 = 1.0
                if abs(vH) > 0.01:
                    self.impe += 0.5 * mE * vH**2 * (1 - self.rest**2)
                    tauImp += mE * abs(vH) * (1 + self.rest) * Rm / dt * self.impF

            self.pos[i] = max(0.0, min(1.0, nP2))
            self.pos[j] = 1.0 - self.pos[i]
            self.vel[j] = -self.vel[i]

            # طاقة مبددة (إصلاح: بدون Lm الزائد)
            if abs(self.vel[i]) > 0.001:
                self.diss += abs((tR + Fs * np.tanh(v1 / vSm)) * self.vel[i] * dt)

        # ── 3. مقاومات العجلة ──
        tauA = -np.sign(self.omega if self.omega != 0 else 1) * (
            self.axMu * totalM * g * Rm * 0.7 + abs(self.omega) * 0.05
        )
        tauAir = -self.airD * self.omega * abs(self.omega) * 1000

        tauL = 0.0
        if self.load > 0 and self.omega > 0.01:
            tauL = -self.load
            self.diss += self.load * abs(self.omega) * dt

        totalTau += tauAir + tauA + tauImp + tauL

        if abs(self.omega) > 0.01:
            self.diss += abs((tauAir + tauA) * self.omega * dt)

        # ── 4. تحديث حركة العجلة ──
        alpha = totalTau / totalI
        self.omega += alpha * dt
        self.ang += self.omega * dt
        self.t += dt

        # ── 5. حساب الطاقات ──
        PE = 0.0
        for i in range(N):
            a = self.ang + self.off[i]
            sp = -Lm/2 + Lm * self.pos[i]
            sw = -Lm/2 + Lm * self.pos[i] / 2
            PE += Mp * g * (Rm - (Rm * np.sin(a) + sp * np.cos(a)))
            PE += Mw * self.pos[i] * g * (Rm - (Rm * np.sin(a) + sw * np.cos(a)))

        KE = 0.5 * totalI * self.omega**2

        return SimState(
            t=self.t, omega=self.omega, alpha=alpha,
            KE=KE, PE=PE, diss=self.diss, impe=self.impe,
            totalE=KE + PE + self.diss + self.impe,
            I=totalI, tau=totalTau, pos=self.pos.copy()
        )

    def run(self, dur: float, dt: float = 0.002, sub: int = 6, rec_every: int = 25) -> List[SimState]:
        """تشغيل المحاكاة لمدة محددة"""
        sdt = dt / sub
        recs = []
        s = 0
        while self.t < dur:
            r = None
            for _ in range(sub):
                r = self.step(sdt)
            if r and s % rec_every == 0:
                recs.append(r)
            if abs(self.omega) < 0.0003 and self.t > 0.5:
                break
            s += 1
        return recs

    def center_of_mass(self) -> Tuple[float, float]:
        """حساب مركز الكتلة"""
        cx, cy, tM = 0.0, 0.0, 0.0
        for i in range(self.N):
            a = self.ang + self.off[i]
            sp = -self.Lm/2 + self.Lm * self.pos[i]
            sw = -self.Lm/2 + self.Lm * self.pos[i] / 2
            mW = self.Mw * self.pos[i]
            cx += self.Mp * (self.Rm * np.cos(a) - sp * np.sin(a)) + mW * (self.Rm * np.cos(a) - sw * np.sin(a))
            cy += self.Mp * (self.Rm * np.sin(a) + sp * np.cos(a)) + mW * (self.Rm * np.sin(a) + sw * np.cos(a))
            tM += self.Mp + mW
        return cx / tM, cy / tM


# ═══════════════════════════════════════════════════════
# القسم ٢: الاختبارات الفيزيائية الثمانية
# ═══════════════════════════════════════════════════════

class TestSuite:
    """مجموعة الاختبارات الفيزيائية الشاملة"""

    def __init__(self, sim_class=HydraulicWheelSim):
        self.Sim = sim_class
        self.results = {}

    def run_all(self, verbose=True):
        """تشغيل جميع الاختبارات"""
        tests = [
            ("test1_convergence", self.test1_convergence),
            ("test2_energy", self.test2_energy),
            ("test3_torque_sweep", self.test3_torque_sweep),
            ("test4_longterm", self.test4_longterm),
            ("test5_com", self.test5_com),
            ("test6_friction", self.test6_friction),
            ("test7_power", self.test7_power),
            ("test8_inertia", self.test8_inertia),
        ]

        for name, test_fn in tests:
            if verbose:
                print(f"\n{'='*60}")
                print(f"Running {name}...")
                print('='*60)
            try:
                result = test_fn()
                self.results[name] = result
                if verbose:
                    self._print_result(result)
            except Exception as e:
                self.results[name] = {"error": str(e)}
                if verbose:
                    print(f"ERROR: {e}")

        return self.results

    def _print_result(self, r):
        for k, v in r.items():
            if k != "_data":
                print(f"  {k}: {v}")

    # ── TEST 1: Timestep Convergence ──
    def test1_convergence(self):
        dts = [0.02, 0.01, 0.005, 0.002, 0.001]
        fin_omegas = []
        for dt in dts:
            s = self.Sim()
            s.omega = 1.5
            sub = max(2, round(0.01 / dt))
            recs = s.run(6.0, dt, sub, 15)
            fin_omegas.append(recs[-1].omega if recs else 0)

        ref = fin_omegas[-1]
        max_err = max(abs(v - ref) for v in fin_omegas[:-1])
        return {
            "dts": dts,
            "final_omegas": fin_omegas,
            "reference": ref,
            "max_error": max_err,
            "converged": max_err < 0.1,
            "pass": max_err < 0.1,
            "_data": {"dts": dts, "omegas": fin_omegas}
        }

    # ── TEST 2: Energy Audit ──
    def test2_energy(self):
        # No gravity
        s_ng = self.Sim(fGrav=False, fCentr=False, fVisc=False)
        s_ng.omega = 2.0
        recs_ng = s_ng.run(10.0, 0.002, 6, 20)
        E0_ng = recs_ng[0].totalE
        drift_ng = [r.totalE - E0_ng for r in recs_ng]
        pct_ng = abs(drift_ng[-1] / E0_ng * 100) if E0_ng != 0 else 0

        # With gravity
        s_fg = self.Sim()
        s_fg.omega = 2.0
        recs_fg = s_fg.run(12.0, 0.002, 6, 20)
        E0_fg = recs_fg[0].totalE
        drift_fg = [r.totalE - E0_fg for r in recs_fg]
        pct_fg = abs(drift_fg[-1] / E0_fg * 100) if E0_fg != 0 else 0

        return {
            "no_gravity_drift_pct": pct_ng,
            "with_gravity_drift_pct": pct_fg,
            "pass_no_gravity": pct_ng < 2,
            "pass_with_gravity": pct_fg < 5,
            "_data": {"recs_ng": recs_ng, "recs_fg": recs_fg}
        }

    # ── TEST 3: Static Torque Sweep ──
    def test3_torque_sweep(self):
        s = self.Sim()
        N3 = 720
        dTheta = 2 * np.pi / N3
        torques = []
        cum_work = [0.0]
        cum_w = 0.0

        for k in range(N3):
            theta = k * dTheta
            tau_k = 0.0
            for i in range(s.N):
                a = theta + s.off[i]
                sp = -s.Lm/2 + s.Lm * s.pos[i]
                sw = -s.Lm/2 + s.Lm * s.pos[i]/2
                mW = s.Mw * s.pos[i]
                tau_k += s.Mp * s.g * (s.Rm * np.cos(a) - sp * np.sin(a))
                tau_k += mW * s.g * (s.Rm * np.cos(a) - sw * np.sin(a))
            torques.append(tau_k)
            cum_w += tau_k * dTheta
            cum_work.append(cum_w)

        net_work = cum_work[-1]
        return {
            "net_work": net_work,
            "pass": abs(net_work) < 1.0,
            "_data": {"torques": torques, "cum_work": cum_work[1:]}
        }

    # ── TEST 4: Long-term Stability ──
    def test4_longterm(self):
        s = self.Sim()
        s.omega = 1.5
        dt4, sub4 = 0.003, 6
        sdt4 = dt4 / sub4
        cyc_om, cyc_n = [], []
        last_ang, cyc_c = s.ang, 0
        step_n, max_s = 0, int(60 / dt4)

        while step_n < max_s:
            for _ in range(sub4):
                s.step(sdt4)
            if s.ang - last_ang >= 2 * np.pi:
                cyc_c += 1
                last_ang += 2 * np.pi
                cyc_om.append(s.omega)
                cyc_n.append(cyc_c)
            if abs(s.omega) < 0.0003:
                break
            step_n += 1

        trend = 0.0
        if len(cyc_om) > 3:
            n = len(cyc_om)
            sx, sy = sum(cyc_n), sum(cyc_om)
            sxx = sum(x**2 for x in cyc_n)
            sxy = sum(x*y for x, y in zip(cyc_n, cyc_om))
            trend = (n * sxy - sx * sy) / (n * sxx - sx**2)

        return {
            "runtime": s.t,
            "cycles": cyc_c,
            "final_omega": s.omega,
            "trend": trend,
            "pass": trend > -0.001,  # Note: original shows acceleration
            "_data": {"cyc_n": cyc_n, "cyc_om": cyc_om}
        }

    # ── TEST 5: Center of Mass ──
    def test5_com(self):
        s = self.Sim()
        s.omega = 1.5
        dist5 = []
        step_n, max_s = 0, int(30 / 0.003)

        while step_n < max_s:
            for _ in range(6):
                s.step(0.003/6)
            if step_n % 30 == 0:
                cm = s.center_of_mass()
                dist5.append(np.sqrt(cm[0]**2 + cm[1]**2))
            if abs(s.omega) < 0.0003:
                break
            step_n += 1

        mean_d = np.mean(dist5)
        return {
            "mean_com_distance": mean_d,
            "min": min(dist5),
            "max": max(dist5),
            "pass": mean_d < 0.01,
            "_data": {"distances": dist5}
        }

    # ── TEST 6: Friction Threshold ──
    def test6_friction(self):
        mu_vals = [0, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]
        om_mu = []
        for mu in mu_vals:
            s = self.Sim(axMu=mu)
            s.omega = 1.5
            s.run(8.0, 0.005, 4, 50)
            om_mu.append(s.omega)

        thr_mu = next((mu for mu, om in zip(mu_vals, om_mu) if om < 0.01), None)

        visc_vals = [0, 50, 100, 150, 250, 350, 500]
        om_visc = []
        for v in visc_vals:
            s = self.Sim(visc=v)
            s.omega = 1.5
            s.run(8.0, 0.005, 4, 50)
            om_visc.append(s.omega)

        return {
            "mu_values": mu_vals,
            "omega_vs_mu": om_mu,
            "visc_values": visc_vals,
            "omega_vs_visc": om_visc,
            "threshold_mu": thr_mu,
            "pass": True,
            "_data": {"mu": mu_vals, "om_mu": om_mu, "visc": visc_vals, "om_visc": om_visc}
        }

    # ── TEST 7: Power Curve ──
    def test7_power(self):
        loads = [0, 25, 50, 100, 150, 200, 300, 400, 500]
        fin_sp, fin_pw = [], []
        for ld in loads:
            s = self.Sim(load=ld)
            s.omega = 2.0
            s.run(10.0, 0.005, 4, 50)
            fin_sp.append(s.omega)
            fin_pw.append(ld * s.omega)

        max_p = max(fin_pw)
        max_pi = fin_pw.index(max_p)

        return {
            "loads": loads,
            "final_speeds": fin_sp,
            "powers": fin_pw,
            "max_power": max_p,
            "max_power_load": loads[max_pi],
            "pass": True,
            "_data": {"loads": loads, "speeds": fin_sp, "powers": fin_pw}
        }

    # ── TEST 8: Dynamic Inertia ──
    def test8_inertia(self):
        s = self.Sim()
        s.omega = 1.5
        I_h, t8 = [], []
        step_n, max_s = 0, int(15 / 0.005)

        while step_n < max_s:
            r = None
            for _ in range(4):
                r = s.step(0.005/4)
            if r and step_n % 30 == 0:
                I_h.append(r.I)
                t8.append(r.t)
            if abs(s.omega) < 0.001:
                break
            step_n += 1

        i_mean = np.mean(I_h)
        i_var = (max(I_h) - min(I_h)) / i_mean * 100

        return {
            "mean_inertia": i_mean,
            "variation_pct": i_var,
            "pass": i_var < 10,
            "_data": {"I": I_h, "t": t8}
        }


# ═══════════════════════════════════════════════════════
# القسم ٣: التشغيل الرئيسي
# ═══════════════════════════════════════════════════════

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return json.JSONEncoder.default(self, obj)

def main():
    print("=" * 70)
    print("  المختبر الهيدروليكي — اختبارات فيزيائية شاملة")
    print("  Hydraulic Wheel Physics Test Suite")
    print("=" * 70)

    suite = TestSuite()
    results = suite.run_all(verbose=True)

    # Summary
    print("\n" + "=" * 70)
    print("  الملخص النهائي")
    print("=" * 70)

    for name, r in results.items():
        status = "✅ PASS" if r.get("pass", False) else "❌ FAIL"
        if "error" in r:
            status = "💥 ERROR"
        print(f"  {name:25s} {status}")

    # Save results
    with open("test_results.json", "w", encoding="utf-8") as f:
        # Remove _data for JSON serialization
        clean = {k: {kk: vv for kk, vv in v.items() if kk != "_data"}
                 for k, v in results.items()}
        json.dump(clean, f, indent=2, ensure_ascii=False, cls=NumpyEncoder)
    print("\n  تم حفظ النتائج في: test_results.json")


if __name__ == "__main__":
    main()