# ⚙️ Asymmetric Hydraulic Pendulum Wheel (AHPW)
### عجلة البندول الهيدروليكي غير المتماثل

[![License: AHPW-OIL](https://img.shields.io/badge/License-AHPW--Open--Innovation-blue.svg)](./LICENSE)
[![DOI](https://img.shields.io/badge/DOI-pending%20Zenodo-orange.svg)]()
[![Status](https://img.shields.io/badge/Status-Open%20for%20Experimental%20Verification-green.svg)]()

---

## 🔬 What Is This?

The **Asymmetric Hydraulic Pendulum Wheel (AHPW)** is a novel mechanical concept invented by **Mustafa Sayed Sharif** (Co-Founder, MaelshPro Technical Solutions) in 2026.

The core idea: fluid-coupled piston pairs arranged radially on a rotating wheel continuously transfer fluid mass toward the gravitationally descending (right) side of the wheel — creating a **persistent net torque imbalance** throughout the full 360° rotation cycle.

> **This is not a perpetual motion claim.**
> It is a mechanically coherent concept with positive net gravitational work (∫τ dθ = +4,607 J per revolution in pure-physics simulation) that requires physical experimental verification.

---

## 💡 The Core Innovation

```
Classical overbalanced wheel:   ∫τ dθ = 0  (always fails)
AHPW:                           ∫τ dθ > 0  (simulation result — needs physical test)
```

Three features make AHPW different from all prior attempts:

| Feature | Classical Wheel | AHPW |
|---|---|---|
| Internal state change | None | Continuous fluid redistribution |
| Dead points | Exist | Eliminated (N ≥ 8 chambers) |
| Transfer driver | None | Gravity on piston (M_piston > M_fluid) |
| ∫τ dθ per revolution | = 0 | ≠ 0 (numerically confirmed) |

---

## ⚙️ How It Works

Each **Z-Piston Unit** consists of:
- Two cylindrical chambers at diametrically opposite points on the wheel rim
- A Z-shaped fluid connector linking them at an angle
- Two heavy pistons (80 kg each) sliding inside the chambers
- Fluid (≤ 60 kg) trapped between the pistons

**Key principle:** Since `M_piston (80 kg) >> M_fluid (60 kg)`, gravity drives the pistons to push fluid toward the upper-right chamber — without any external pump.

```
Phase 0°  → Fluid lower-right  → Right side heavier → Wheel tilts right
Phase 45° → Pistons migrate right, fluid rises through Z-channel
Phase 90° → Fluid upper-right  → Peak torque: ~647 N·m
Phase 135°→ Fluid held right   → Continuous positive torque
Phase 180°→ Next unit activates → No dead point, no interruption
```

With **8 chambers** offset at 45° intervals, at least one unit is always actively pushing fluid rightward. The right side of the wheel is **always heavier**.

---

## 📊 Simulation Results

Simulation parameters:
- N = 8 chambers, M_p = 80 kg, M_w_max = 60 kg, R = 0.30 m
- Fluid viscosity η = 150 N·s/m (intentionally high for stress-test)
- Impact factor = **0.00** (pure physics — no artificial energy injection)

| Quantity | Value |
|---|---|
| ω start | 0.500 rad/s |
| ω after 1 revolution | 1.116 rad/s ✅ |
| ∫τ dθ — gravitational work | **+4,607 J** |
| Net mean torque | 733 N·m |
| Axle friction torque | ~10.2 N·m |
| Safety margin | **72×** |

> The wheel **accelerates** under pure gravitational physics with no artificial energy injection.

---

## 🔧 Build Your Own — Minimum Viable Prototype

| Specification | Value |
|---|---|
| Number of chambers | 8 |
| Orbital radius R | 0.30 – 0.60 m |
| Piston mass | ≥ 5× fluid mass |
| Fluid | Distilled water (η ≈ 0.001 N·s/m) |
| Chamber material | PVC or aluminum tube |
| Axle bearing | Ball bearings (μ < 0.005) |
| Z-connector angle | 30° – 45° from radial axis |

**Test protocol:**
1. Apply initial push: ω₀ = 0.5 rad/s
2. Record ω(t) for 120+ seconds
3. If dω/dt > 0 → positive result
4. If wheel maintains speed with load → partial positive result
5. If wheel stops → optimize M_p/M_w ratio and reduce friction

---

## 📁 Repository Contents

```
/
├── README.md                          ← This file
├── LICENSE                            ← Open Innovation License
├── papers/
│   ├── AHPW_Scientific_Paper_EN.docx ← Full English scientific paper
│   ├── AHPW_Scientific_Paper_AR.docx ← Full Arabic scientific paper
│   ├── AHPW_Detailed_Academic_EN.docx← Detailed academic paper (EN)
│   ├── AHPW_Detailed_Academic_AR.docx← Detailed academic paper (AR)
│   └── AHPW_Simple_Guide_AR.docx     ← Simplified guide for builders
├── simulation/
│   └── simulation-v15.html           ← Interactive browser simulation
├── results/
│   ├── energy_profile.png            ← Energy balance charts
│   └── one_cycle_analysis.png        ← Single revolution torque analysis
└── diagrams/
    └── (phase diagrams — coming soon)
```

---

## 🌍 Open Questions for the Scientific Community

- Does ∫τ dθ > 0 persist across many consecutive revolutions?
- What is ω_max before fluid transfer lag eliminates the asymmetry?
- What is the optimal M_p/M_w ratio as a function of angular velocity?
- Can the mechanism sustain useful mechanical output power?
- How do results change with water-accurate viscosity (η = 0.001)?

**If you build a prototype — please document and share your results openly.**

---

## 👤 Inventor

**Mustafa Sayed Sharif**
Co-Founder, MaelshPro Technical Solutions
📧 maelshspro@gmail.com

---

## 📜 License

This project is licensed under the **AHPW Open Innovation License**.
See [LICENSE](./LICENSE) for full terms.

**Short version:**
- ✅ Free to study, build, test, and share for peaceful purposes
- ✅ Academic and educational use — no permission needed
- ❌ Military or harmful applications — strictly prohibited
- ❌ Commercial exploitation without written consent from inventor

---

## 🤝 Contributing

Contributions welcome:
- Physical prototype builds and test results
- Improved simulation models
- Mathematical analysis and peer review
- Translations

Please open an Issue or Pull Request.
If you build a prototype, open an Issue titled **"Prototype Build — [Your Country]"** and share your measurements.

---

*"The question is not whether perpetual motion is achievable — it is whether a sustained torque differential can be maintained long enough to be practically useful."*
*— Mustafa Sayed Sharif, 2026*
