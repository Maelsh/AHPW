import numpy as np
import matplotlib.pyplot as plt

# إعدادات المحاكاة
n_chambers = 12
m_piston = 80    # kg
m_water_max = 60 # kg
g = 9.81
radius = 1.0     # meter
dt = 0.01
steps = 3000

# المصفوفات لتخزين النتائج
angles = [0.0]
omega = [0.0]
torques = []

# المحاكاة الفيزيائية
for i in range(steps):
    theta = angles[-1]
    current_torque = 0
    inertia = 500 # زيادة القصور الذاتي للهيكل

    for j in range(n_chambers):
        # الزاوية الخاصة بكل حجرة
        angle_j = theta + (j * 2 * np.pi / n_chambers)

        # محاكاة هيدروليكية: الحجرات بين زاوية 0 و 180 درجة (الجهة اليمنى) تحتوي على ماء
        # الحجرات في الجهة الأخرى فارغة
        effective_angle = angle_j % (2 * np.pi)
        if 0 <= effective_angle <= np.pi:
            current_mass = m_piston + m_water_max
        else:
            current_mass = m_piston

        # العزم = القوة * الذراع الأفقي (sin theta)
        current_torque += current_mass * g * np.sin(effective_angle) * radius
        inertia += current_mass * radius**2

    # إضافة مقاومة بسيطة (friction)
    current_torque -= 0.5 * omega[-1]

    # حساب التسارع الزاوي
    alpha = current_torque / inertia

    # تكامل يويلر
    new_omega = omega[-1] + alpha * dt
    new_angle = angles[-1] + new_omega * dt

    omega.append(new_omega)
    angles.append(new_angle)
    torques.append(current_torque)

# تحليل النتائج
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

ax1.plot(np.arange(len(omega))*dt, omega, color='blue')
ax1.set_title("تطور السرعة الزاوية مع توزيع المياه المتغير")
ax1.set_ylabel("rad/s")
ax1.set_xlabel("الزمن (ثانية)")
ax1.grid(True)

ax2.plot(np.arange(len(torques))*dt, torques, color='red')
ax2.set_title("العزم الكلي المتولد (Net Torque)")
ax2.set_xlabel("الزمن (ثانية)")
ax2.set_ylabel("N.m")
ax2.grid(True)

plt.tight_layout()
plt.show()

print(f"أقصى سرعة دورانية وصل إليها النظام: {max(omega):.4f} rad/s")
print(f"العزم النهائي: {torques[-1]:.2f} N.m")