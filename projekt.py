import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

# Domyślne parametry układu
ampl = 1.0 #częstotliwość
wyp = 0.5 #wypełnienie
f = 10

M = 1.0  # masa wózka
k = 10.0  # stała sprężyny
b = 1.0  # współczynnik tłumienia
dt = 0.001   # krok czasowy, 
T = 20.0  # czas symulacji

N = int(T / dt)  # liczba kroków symulacji
force_type = "step"
t0 = 1.0          # początek impulsu [s]
T_pulse = 2.0     # długość impulsu [s]
def force_signal(t, ampl, wyp, f, type="step"):  #SYGNAŁ WEJŚCIOWY
    okres  = 1.0/f
    tt=t % okres
    if type == "step":
        return ampl if t > 1 else 0.0
    elif type == "triangle":
        if tt < wyp * okres:
            return ampl * (tt / (wyp * okres))  # Narastanie
        else:
            return ampl * (1 - (tt - wyp * okres) / ((1 - wyp) * okres))  # Opadanie
    elif type == "sine":
        return ampl*np.sin(2 * np.pi * f * t)
    elif type == "square":
        return ampl if  tt <= okres*wyp else 0.0
    elif type =="pulse":
        return ampl if (t >= t0) and (t < t0 + T_pulse) else 0.0
    return 0

def euler_method(N, dt, M, k, b, ampl, wyp, f, force_type):
    x, v = 0.0, 0.0  # Resetowanie wartości początkowych
    positions, velocities = [], []
    for i in range(N):
        t = i * dt
        F = force_signal(t, ampl, wyp, f, force_type)
        a = (F - b * v - k * x) / M
        v += a * dt
        x += v * dt
        positions.append(x)
        velocities.append(v)
    return np.array(positions), np.array(velocities)

def rk4_method(N, dt, M, k, b, ampl, wyp, f, force_type):
    x, v = 0.0, 0.0  # Resetowanie wartości początkowych
    positions, velocities = [], []
    for i in range(N):
        t = i * dt


        def acceleration(x, v, t):
            F = force_signal(t, ampl, wyp, f, force_type)
            return (F - b * v - k * x) / M

        k1v = acceleration(x, v, t) * dt
        k1x = v * dt

        k2v = acceleration(x + k1x / 2, v + k1v / 2, t + dt / 2) * dt
        k2x = (v + k1v / 2) * dt

        k3v = acceleration(x + k2x / 2, v + k2v / 2, t + dt / 2) * dt
        k3x = (v + k2v / 2) * dt

        k4v = acceleration(x + k3x, v + k3v, t + dt) * dt
        k4x = (v + k3v) * dt

        v += (k1v + 2 * k2v + 2 * k3v + k4v) / 6
        x += (k1x + 2 * k2x + 2 * k3x + k4x) / 6

        positions.append(x)
        velocities.append(v)
    return np.array(positions), np.array(velocities)

def plot_simulation():
    global M, k, b, force_type
    x_euler, v_euler = euler_method(N, dt, M, k, b, ampl, wyp, f, force_type)
    x_rk4, v_rk4 = rk4_method(N, dt, M, k, b, ampl, wyp, f, force_type)
    times = np.linspace(0, T, N)
    sila_wej = np.array([force_signal(t, ampl, wyp, f, force_type) for t in times])

    ax1 = plt.subplot(3, 1, 1)  # Położenie
    ax2 = plt.subplot(3, 1, 2)  # Prędkość
    ax3 = plt.subplot(3, 1, 3)  # Siła

    ax1.clear()
    ax2.clear()
    ax3.clear()

    ax1.plot(times, x_euler, label="Euler")
    ax1.plot(times, x_rk4, label="RK4", linestyle='dashed')
    ax1.set_title("Położenie wózka")
    ax1.set_ylabel("Położenie x [m]")
    ax1.legend()
    ax1.grid(True)

    ax2.plot(times, v_euler, label="Euler")
    ax2.plot(times, v_rk4, label="RK4", linestyle='dashed')
    ax2.set_title("Prędkość wózka")
    ax2.set_ylabel("Prędkość v [m/s]")
    ax2.legend()
    ax2.grid(True)

    ax3.plot(times, sila_wej, label=f"Siła ({force_type})", color='red')
    ax3.set_title("Siła wejściowa")
    ax3.set_xlabel("Czas [s]")
    ax3.set_ylabel("Siła F [N]")
    ax3.legend()
    ax3.grid(True)

    plt.draw()

def update(val):
    global M, k, b,f,ampl,wyp,t0,T_pulse,dt
    M = slider_M.val
    k = slider_k.val
    b = slider_b.val
    f = slider_f.val
    wyp = slider_wyp.val
    ampl = slider_ampl.val
    t0 = slider_t0.val
    T_pulse = slider_pulse.val
    plot_simulation()

def update_force(label):
    global force_type
    force_type = label
    plot_simulation()

fig = plt.figure(figsize=(11, 11))
plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.3, hspace=0.28)

slider_height = 0.03
slider_start = 0.21
slider_spacing = 0.027

ax_M = plt.axes([0.15, slider_start, 0.35, slider_height])
ax_k = plt.axes([0.15, slider_start - slider_spacing, 0.35, slider_height])
ax_b = plt.axes([0.15, slider_start - slider_spacing*2, 0.35, slider_height])
ax_f = plt.axes([0.15, slider_start - slider_spacing*3, 0.35, slider_height])
ax_ampl = plt.axes([0.15, slider_start - slider_spacing*4, 0.35, slider_height])
ax_wyp = plt.axes([0.15, slider_start - slider_spacing*5, 0.35, slider_height])
ax_t0 = plt.axes([0.15, slider_start - slider_spacing*6, 0.35, slider_height])
ax_pulse = plt.axes([0.15, slider_start - slider_spacing*7, 0.35, slider_height])

radio_ax = plt.axes([0.6, 0.15 - slider_spacing*2, 0.25, slider_height*4])

slider_M = Slider(ax_M, 'Masa (M)', 0.1, 5.0, valinit=M)
slider_k = Slider(ax_k, 'Sprężystość (k)', 1.0, 50.0, valinit=k)
slider_b = Slider(ax_b, 'Tłumienie (b)', 0.1, 5.0, valinit=b)
slider_f = Slider(ax_f, 'Częstotliwość', 0.1, 20.0, valinit=f)
slider_ampl = Slider(ax_ampl, 'Amplituda', 0.1, 10.0, valinit=ampl)
slider_wyp = Slider(ax_wyp, 'Wypełnienie', 0.01, 0.99, valinit=wyp)
slider_t0    = Slider(ax_t0,    'Start impulsu t₀',    0.0, T,    valinit=t0)
slider_pulse = Slider(ax_pulse, 'Czas trwania Tₚ',     0.0, T,    valinit=T_pulse)

radio_force = RadioButtons(radio_ax, ('step','pulse','triangle','sine','square'))


slider_M.on_changed(update)
slider_k.on_changed(update)
slider_b.on_changed(update)
slider_f.on_changed(update)
slider_ampl.on_changed(update)
slider_wyp.on_changed(update)
slider_t0.on_changed(update)
slider_pulse.on_changed(update)
radio_force.on_clicked(update_force)


plot_simulation()
plt.show()
