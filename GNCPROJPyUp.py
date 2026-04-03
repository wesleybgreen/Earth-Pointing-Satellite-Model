import numpy as np
import scipy.linalg
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from tqdm import tqdm

# GNC Earth-Pointing project

# Initial Variables and Constants
quatin = np.sqrt(2)/2 * np.array([1, 0, 0, 1])
quat0 = quatin/np.linalg.norm(quatin)
w0 = np.array([13.53*np.pi/180, 13.53*np.pi/180, 13.053*np.pi/180])
hrw0 = np.array([0, 0, 0])
J = np.array([[0.035, -0.001, 0.0004], [-0.001, 0.035, -0.0008], [0.0004, -0.0008, 0.006]])
Jinv = np.linalg.inv(J)

pi = np.pi
sqrt = np.sqrt
kep0 = np.array([15.49442598, 0.0007785, 51.6334*pi/180, 312.1983*pi/180, 38.3265*pi/180, 321.8276*pi/180, .00016677])
ain = ((3.986e5)/((kep0[0]*2*pi/24/60/60)**2))**(1/3)
T = 2*pi*sqrt((ain**3)/(3.986e5))

# initial state
s0 = np.concatenate([quat0.reshape(-1, 1), w0.reshape(-1, 1), hrw0.reshape(-1, 1)])
tspan = [0, T]
teval = np.linspace(0, T, 20500)



def dynamics(t, sflat, J, kep0):
    
    s = sflat.reshape(-1,1)
    rng = np.random.default_rng()

    # Extracting state
    quatnew = s[0:4]
    quatcur = quatnew/np.linalg.norm(quatnew)
    wcur = s[4:7]
    hrwcur = s[7:10]

    # Getting current position and vel vector
    rvec, vvec = orbitprop(t, kep0)

    # Calculate frame change and error values
    quaterr, werr, RBL = errorfunc(t, rvec, vvec, quatcur, wcur, kep0)

    # Control Law
    tauctrl = ctrllaw(quaterr, werr)

    # Reaction Wheels
    BECI = magfield(t, rvec).reshape(3,1)
    BN = BECI/np.linalg.norm(BECI)
    BBtrue = (quat2dcm(quatcur) @ BN).reshape(3,1)
    BBmeas = BBtrue + rng.normal(0, 1.5*pi/180, (3,1))
    BBmeas /= np.linalg.norm(BBmeas)
    tauprod, hrwdot, md = RW(-tauctrl, wcur, hrwcur, BBtrue)

    # External Torque
    tauext = exttorque(rvec, vvec, J, RBL)

    wdot = (Jinv) @ (tauext.reshape(-1,1) -tauprod - np.cross(wcur, J @ wcur, axis=0))
    quatdot = kinematics(quatcur, wcur)

    sdot = np.concatenate([quatdot.reshape(-1,1), wdot.reshape(-1,1), hrwdot.reshape(-1,1)])
    return sdot.flatten()

def orbitprop(t, kep):
    mu = 3.986e5
    n = kep[0]*2*pi/24/60/60
    e = kep[1]
    i = kep[2]
    RAAN = kep[3]
    AOP = kep[4]
    Min = kep[5]
    ndot = kep[6]*2*pi/(86400**2)

    a = (mu/n**2)**(1/3)
    M = np.mod(Min + n*t + ndot*t**2, 2*pi)
    E = M

    for j in range(1,10):
        fe = E - e*np.sin(E) - M
        df = 1 - e*np.cos(E)
        E = E - fe/df
        if np.max(np.abs(fe)) < 1e-12:
            break

    f = 2 * np.arctan2(sqrt(1+e)*np.sin(E/2), sqrt(1-e)*np.cos(E/2))
    r = (a*(1-e**2))/(1+e*np.cos(f))
    sw, cw = np.sin(AOP), np.cos(AOP)
    sr, cr = np.sin(RAAN), np.cos(RAAN)
    si, ci = np.sin(i), np.cos(i)
    R_ECI = np.array([[cw*cr-sr*sw*ci, -cr*sw - sr*cw*ci, sr*si],
                      [sr*cw + cr*sw*ci, -sr*sw + cr*cw*ci, -cr*si], 
                      [sw*si, cw*si, ci]])

    r_vec_peri = np.array([r*np.cos(f), r*np.sin(f), 0])
    v_vec_peri = np.array([-sqrt(mu/(a*(1-e**2)))*np.sin(f), sqrt(mu/(a*(1-e**2)))*(e+np.cos(f)), 0])

    r_vec_ECI = R_ECI @ r_vec_peri.reshape(-1,1)
    v_vec_ECI = R_ECI @ v_vec_peri.reshape(-1,1)

    r_vec = r_vec_ECI.reshape(-1,1)
    v_vec = v_vec_ECI.reshape(-1,1)

    return r_vec.reshape(-1,1), v_vec.reshape(-1,1)

def errorfunc(t, rvec, vvec, quatBI, wBI, kep):
    zLVLH = -rvec/np.linalg.norm(rvec)
    yLVLH = -np.cross(rvec, vvec, axis=0)/np.linalg.norm(np.cross(rvec, vvec, axis=0))
    xLVLH = np.cross(yLVLH, zLVLH, axis=0)

    RLI = np.column_stack((xLVLH, yLVLH, zLVLH))
    RBI = (quat2dcm(quatBI))
    RBL = RBI @ RLI
    quaterr = (dcm2quat(RBL))
    
    #if quaterr[0] < 0:
    #    quaterr = -quaterr

    ain = ((3.986e5)/((kep[0]*2*pi/24/60/60)**2))**(1/3)
    wLI = np.array([0, -sqrt(3.986e5/(ain**3)), 0])
    werr = wBI.reshape(-1,1) - RBL @ wLI.reshape(-1,1)

    return quaterr.reshape(-1,1), werr.reshape(-1,1), RBL

def ctrllaw(quaterr, werr):
    Kp = 0.005
    Kd = 0.005

    if quaterr[0] < 0:
        quaterr = -quaterr

    tauctrl = -Kp*quaterr[1:4] - Kd*werr
    return tauctrl.reshape(-1,1)

def RW(taucom, wcur, hrwcur, BB):
    taucom = taucom.reshape(-1,1)
    wcur = wcur.reshape(-1,1)
    hrwcur = hrwcur.reshape(-1,1)
    BB = BB.reshape(-1,1)

    BBmag = np.linalg.norm(BB)
    if BBmag > 1e-12:
        mdcom = -0.01 * np.cross(BB, hrwcur, axis=0) / BBmag**2
        md = np.clip(mdcom, -0.2, 0.2)
    else:
        md = np.zeros((3,1))

    taumag = np.cross(md, BB, axis=0)
    taurw = taucom - taumag

    taucomsat = np.clip(taurw, -0.002, 0.002)

    hrwdot = taucomsat

    for i in range(3):
        if hrwcur[i] >= 0.015 and hrwdot[i] > 0:
            hrwdot[i] = 0
        elif hrwcur[i] <= -0.015 and hrwdot[i] < 0:
            hrwdot[i] = 0
    
    tauprod = hrwdot.reshape(-1,1) + np.cross(wcur.reshape(-1,1), hrwcur.reshape(-1,1), axis=0)

    return tauprod.reshape(-1,1), hrwdot.reshape(-1,1), md
            
def exttorque(rvecECI, vvecECI, J, RBL):
    mu = 3.986e5

    # Gravity Gradient
    rvec = RBL @ rvecECI.reshape(-1,1)
    rmag = np.linalg.norm(rvec)
    n = -rvec/rmag

    taugg = np.cross(3*mu/(rmag**3)*n, J @ n, axis=0)

    # Aerodynamic Torque
    we = 0.000072921158553*np.array([0,0,1]).reshape(-1,1)
    vrelI = vvecECI.reshape(-1,1) + np.cross(we, rvecECI.reshape(-1,1), axis=0)
    vrelB = RBL @ vrelI

    costh = np.dot(n.T, vrelB)/ np.linalg.norm(vrelB)

    h = np.linalg.norm(rvec) - 6371
    p = density(h)
    Si = 2*0.01 + 4*0.03
    ri = np.array([0.05, 0, -0.01]).reshape(-1,1)

    F = -0.5*p*2.2*np.linalg.norm(vrelB)*vrelB*max(costh, 0)*Si
    tauaero = np.cross(ri, F, axis=0)

    return taugg.reshape(-1,1) + tauaero.reshape(-1,1)

def density(h):
    if h >= 0 and h <= 25:
        h0, po, H = 0, 1.225, 8.44
    elif h > 25 and h <= 30:
        h0, po, H = 25, 3.899e-2, 6.49
    elif h > 30 and h <= 35:
        h0, po, H = 30, 1.774e-2, 6.75
    elif h > 35 and h<= 40:
        h0, po, H = 35, 8.279e-3, 7.07
    elif h > 40 and h <= 45:
        h0, po, H = 40, 3.972e-3, 7.47
    elif h > 45 and h <= 50:
        h0, po, H = 45, 1.995e-3, 7.83
    elif h > 50 and h <= 55:
        h0, po, H =50, 1.057e-3, 7.95
    elif h > 55 and h <= 60:
        h0, po, H =55, 5.821e-4, 7.73
    elif h > 60 and h <= 65:
        h0, po, H =60, 3.206e-4, 7.29
    elif h > 65 and h <= 70:
        h0, po, H =65, 1.718e-4, 6.81
    elif h > 70 and h <= 75:
        h0, po, H =70, 8.770e-5, 6.33
    elif h > 75 and h <= 80:
        h0, po, H =75, 4.178e-5, 6
    elif h > 80 and h <= 85:
        h0, po, H =80, 1.905e-5, 5.70
    elif h > 85 and h <= 90:
        h0, po, H =85, 8.337e-6, 5.41
    elif h > 90 and h <= 95:
        h0, po, H =90, 3.396e-6, 5.38
    elif h > 95 and h <= 100:
        h0, po, H =95, 1.343e-6, 5.74
    elif h > 100 and h <= 110:
        h0, po, H =100, 5.297e-7, 6.15
    elif h > 110 and h <= 120:
        h0, po, H =110, 9.661e-8, 8.06
    elif h > 120 and h <= 130:
        h0, po, H =120, 2.438e-8, 11.6
    elif h > 130 and h <= 140:
        h0, po, H =130, 8.484e-9, 16.1
    elif h > 140 and h <= 150:
        h0, po, H =140, 3.845e-9, 20.6
    elif h > 150 and h <= 160:
        h0, po, H =150, 2.070e-9, 24.6
    elif h > 160 and h <= 180:
        h0, po, H =160, 1.224e-9, 26.3    
    elif h > 180 and h <= 200:
        h0, po, H =180, 5.464e-10, 33.2
    elif h > 200 and h <= 250:
        h0, po, H =200, 2.789e-10, 38.5
    elif h > 250 and h <= 300:
        h0, po, H =250, 7.248e-11, 46.9
    elif h > 300 and h <= 350:
        h0, po, H =300, 2.418e-11, 52.5
    elif h > 350 and h <= 400:
        h0, po, H =350, 9.158e-12, 56.4
    elif h > 400 and h <= 450:
        h0, po, H =400, 3.725e-12, 59.4
    elif h > 450 and h <= 500:
        h0, po, H =450, 1.585e-12, 62.2
    elif h > 500 and h <= 600:
        h0, po, H =500, 6.967e-13, 65.8
    elif h > 600 and h <= 700:
        h0, po, H =600, 1.454e-13, 79
    elif h > 700 and h <= 800:
        h0, po, H =700, 3.614e-14, 109
    elif h > 800 and h <= 900:
        h0, po, H =800, 1.170e-14, 164
    elif h > 900 and h <= 1000:
        h0, po, H =900, 5.245e-15, 225
    elif h > 1000:
            h0, po, H =1000, 3.019e-15, 268

    return po*np.exp(-(h-h0)/H)

def kinematics(quat, w):
    Omega = np.array([[0, -w[0,0], -w[1,0], -w[2,0]],
                     [w[0,0], 0, w[2,0], -w[1,0]],
                     [w[1,0], -w[2,0], 0, w[0,0]],
                     [w[2,0], w[1,0], -w[0,0], 0]])
    dq = 0.5 * Omega @ quat
    return dq.reshape(-1,1)

def quat2dcm(quat):
    quat = quat/np.linalg.norm(quat)
    w, x, y, z = quat.flatten()

    # Precompute terms (x2 is 2*x, etc.)
    x2, y2, z2 = x*2, y*2, z*2
    xx, xy, xz = x*x2, x*y2, x*z2
    yy, yz, zz = y*y2, y*z2, z*z2
    wx, wy, wz = w*x2, w*y2, w*z2
    
    # Standard Inertial-to-Body DCM
    dcm = np.array([
        [1.0 - (yy + zz), xy + wz, xz - wy],
        [xy - wz, 1.0 - (xx + zz), yz + wx],
        [xz + wy, yz - wx, 1.0 - (xx + yy)]
    ])
    return dcm

def dcm2quat(dcm):
    qo = 0.5*sqrt(np.trace(dcm) + 1)
    q1 = (dcm[1,2] - dcm[2,1])/(4*qo)
    q2 = (dcm[2,0] - dcm[0,2])/(4*qo)
    q3 = (dcm[0,1] - dcm[1,0])/(4*qo)

    return np.array([qo, q1, q2, q3]).reshape(-1,1)

def KFdynamics(t, sflat, J, kep0, tauctrl):
    
    s = sflat.reshape(-1,1)
    rng = np.random.default_rng()

    # Extracting state
    quatnew = s[0:4]
    quatcur = quatnew/np.linalg.norm(quatnew)
    wcur = s[4:7]
    hrwcur = s[7:10]

    # Getting current position and vel vector
    rvec, vvec = orbitprop(t, kep0)

    # Calculate frame change and error values
    quaterr, werr, RBL = errorfunc(t, rvec, vvec, quatcur, wcur, kep0)

    # Control Law
    #tauctrl = ctrllaw(quaterr, werr)

    # Reaction Wheels
    BECI = magfield(t, rvec).reshape(3,1)
    BN = BECI/np.linalg.norm(BECI)
    BBtrue = (quat2dcm(quatcur) @ BN).reshape(3,1)
    BBmeas = BBtrue + rng.normal(0, 1.5*pi/180, (3,1))
    BBmeas /= np.linalg.norm(BBmeas)
    tauprod, hrwdot, md = RW(-tauctrl, wcur, hrwcur, BBtrue)

    # External Torque
    tauext = exttorque(rvec, vvec, J, RBL)

    wdot = (Jinv) @ (tauext.reshape(-1,1) -tauprod - np.cross(wcur, J @ wcur, axis=0))
    quatdot = kinematics(quatcur, wcur)

    sdot = np.concatenate([quatdot.reshape(-1,1), wdot.reshape(-1,1), hrwdot.reshape(-1,1)])
    return sdot.flatten()

def quatmul(q, p):
    q = q.flatten()
    p = p.flatten()
    w1, x1, y1, z1 = q
    w2, x2, y2, z2 = p
    
    res = np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])
    return res.reshape(-1, 1)

def ssvec(t):
    lam0 = 0
    tday = t/(86400)
    ep = 23.44*pi/180

    # sun rotates 360 degrees per year
    lamS = lam0 + (360/365.25)*tday
    lams = lamS*pi/180

    sN = np.array([[np.cos(lams), np.sin(lams)*np.cos(ep), np.sin(lams)*np.sin(ep)]])
    return sN.reshape(3,1)

def skew(a):

    return np.array([[0, -a[2,0], a[1,0]],
    [a[2,0], 0, -a[0,0]],
    [-a[1,0], a[0,0], 0]])

def magfield(t, rvec):
    m = (6371**3)*30034e-9
    thm = 169.7 * pi/180
    alm = 108.2 * pi/180
    mECEF = m * np.array([(np.sin(thm) * np.cos(alm)), (np.sin(thm) * np.sin(alm)), (np.cos(thm))]).reshape(3,1)
    we = 0.000072921158553
    gst = we * t
    R = np.array([[np.cos(gst), -np.sin(gst), 0], [np.sin(gst),  np.cos(gst), 0], [0, 0, 1]])
    mhat = (R @ mECEF).reshape(3,1)
    rmag = np.linalg.norm(rvec)
    B = (1/rmag**5) * (3 * np.dot(mhat.T, rvec)*rvec - rmag**2 * mhat)
    return B.reshape(3,1)

s0flat = s0.flatten()

sol = solve_ivp(
    fun = dynamics, 
    t_span = tspan, 
    y0 = s0flat, 
    t_eval=teval, 
    args=(J, kep0),
    method='RK45'
)

quats = sol.y[0:4, :]
omegas = sol.y[4:7, :]
h_wheels = sol.y[7:10, :]

quat_err_hist = []
w_err_hist = []

for i in range(len(sol.t)):
    r, v = orbitprop(sol.t[i], kep0)
    q_e, w_e, RBL = errorfunc(sol.t[i], r, v, quats[:, i], omegas[:, i], kep0)
    quat_err_hist.append(q_e)
    w_err_hist.append(w_e)

quat_err_hist = np.array(quat_err_hist).squeeze()
w_err_hist = np.array(w_err_hist).squeeze()


time = sol.t

plt.figure()
plt.subplot(3, 1, 1)
plt.plot(time, quat_err_hist)
plt.xlabel("Time (seconds)")
plt.ylabel("Quaternion Error")

plt.subplot(3, 1, 2)
plt.plot(time, w_err_hist)
plt.xlabel("Time (seconds)")
plt.ylabel("Angular Velocity Error")

plt.subplot(3, 1, 3)
plt.plot(time, h_wheels.T)
plt.xlabel("Time (seconds)")
plt.ylabel("Wheel Angular Momentum")

plt.show(block=False)
plt.pause(0.1)


# Set up MEKF

# Initialize variables

tauctrl = np.array([0, 0 , 0])
xtrue = np.zeros((len(sol.t), 10))
xtrue[0, :] = s0flat

quatetrueplt = []
wetrueplt = []
currentstateplt = s0flat

# Initial Error
r0, v0 = orbitprop(sol.t[0], kep0)
q_err0, w_err0, RBL0 = errorfunc(sol.t[0], r0, v0, quat0.reshape(-1,1), w0.reshape(-1,1), kep0)
quatetrueplt = [q_err0.flatten()]
wetrueplt = [w_err0.flatten()]
hrwtrueplt = [hrw0.flatten()]

quateKFplt = []
weKFplt = []
yplt = []

# MEKF Gains

# Gyro Bias
rng = np.random.default_rng()
sigv = 0.001
sigu = 1e-6

etav = rng.normal(0, sigv)
etau = rng.normal(0, sigu)
beta = 0.01 * (pi/180) * (1/3600)
betatrue = np.array([beta, beta, beta]).reshape(-1,1)
beta_est_prior = np.zeros((3, 1))

Q_state = (((6/3600)*(np.pi/180))**2)*np.eye(3)
bias_cov = (1e-2)**2 * np.eye(3)
P_prior = np.block([
    [Q_state, np.zeros((3, 3))],
    [np.zeros((3, 3)), bias_cov]])
H = np.hstack((np.eye(3), np.zeros((3, 3))))
R = np.block([
    [((4.85e-5)**2)*np.eye(3), np.zeros((3,3))],
    [np.zeros((3,3)), ((0.1*pi/180)**2)*np.eye(3)]])
Q = np.block([
    [(sigv**2)*np.eye(3), np.zeros((3, 3))],
    [np.zeros((3, 3)), ((sigu)**2)*np.eye(3)]
])

P = np.zeros((6, 6, len(sol.t)))
P[:, :, 0] = P_prior
betahist = []

xhat = np.zeros((len(sol.t), 7))
q0noisy = quat0.reshape(-1, 1) + 0.0001 * np.random.randn(4, 1)
q0noisy = q0noisy/np.linalg.norm(q0noisy)

xhat[0, :] = np.concatenate([q0noisy, beta_est_prior]).flatten()

for j in tqdm(range(1, len(sol.t)), desc="Simulating"):
    truespan = [sol.t[j-1], sol.t[j]]
    dt = sol.t[j] - sol.t[j-1]

    Q_discrete = np.block([
    [(sigv**2) * dt * np.eye(3), np.zeros((3, 3))],
    [np.zeros((3, 3)), (sigu**2) * dt * np.eye(3)]
    ])
    
    soltrue = solve_ivp(
        fun = KFdynamics,
        t_span = truespan,
        y0 = xtrue[j-1, :],
        t_eval = [sol.t[j]],
        args = (J, kep0, tauctrl),
        method = "RK45"
    )

    xtrue[j, :] = soltrue.y[:, -1]
    quatnew = xtrue[j, 0:4].reshape(-1,1)
    quatcur = quatnew/np.linalg.norm(quatnew)
    wcur = xtrue[j, 4:7].reshape(-1,1)
    hrwcur = xtrue[j, 7:10].reshape(-1,1)

    # Current orbit pos and vel
    rvec, vvec = orbitprop(sol.t[j], kep0)

    # Calculate frame change and error values
    quatetrue, wetrue, RBL = errorfunc(sol.t[j], rvec, vvec, quatcur, wcur, kep0)
    quatetrueplt.append(quatetrue.flatten())
    wetrueplt.append(wetrue.flatten())
    hrwtrueplt.append(hrwcur.flatten())

    # MEKF

    Hlist = []
    Rlist = []
    ylist = []

    # Get the true measurements for current time steps and add noise
    betatrue = betatrue.reshape(3,1) + rng.normal(0, sigu*sqrt(dt), (3,1))
    betahist.append(betatrue.flatten())
    wtrue = xtrue[j, 4:7].reshape(-1,1)
    gyronoise = rng.normal(0, sigv / sqrt(dt), (3, 1))
    wmeas = wtrue + betatrue + gyronoise

    west = wmeas - xhat[j-1, 4:7].reshape(-1,1)

    # Predict new quaternion
    dq = kinematics(xhat[j-1, 0:4].reshape(-1,1), west)
    qpred = xhat[j-1, 0:4].reshape(-1,1) + dq * dt
    qpred /= np.linalg.norm(qpred)
    sigroll = rng.normal(0, (40 * pi/(180*3600)))
    sigbs = rng.normal(0, (5 * pi/(180*3600)))

    dqnoise = np.array([1, 0.5*sigbs, 0.5*sigbs, 0.5*sigroll]).reshape(-1,1)
    dqnoise /= np.linalg.norm(dqnoise)
    qmeas = quatmul(dqnoise, quatcur)

    # Predict Bias
    betapred = xhat[j-1, 4:7].reshape(-1,1)

    # Priori estimate
    xest = np.concatenate([qpred, betapred])

    # Kalman
    
    # Star Tracker
    qestinv = np.array([qpred[0], -qpred[1], -qpred[2], -qpred[3]])
    ydq = quatmul(qestinv, qmeas)
    dthtrue = 2*ydq[1:4].reshape(-1,1)
    #y = np.vstack((dthtrue, ysun))
    #yplt.append(y.flatten())
    Hstar = np.hstack([np.eye(3), np.zeros((3,3))])
    Rstarvec = np.array([(5 * pi/(180*3600))**2, (5 * pi/(180*3600))**2, (40 * pi/(180*3600))**2])
    Rstar = np.diag(Rstarvec)
    #Hlist.append(Hstar)
    #ylist.append(dthtrue)
    #Rlist.append(Rstar)

    # Sun Sensor
    sN = ssvec(sol.t[j]) # Sun pos unit vec
    Aqpred = quat2dcm(qpred).reshape(3,3)
    Aqtrue = quat2dcm(quatcur).reshape(3,3)
    sBest = (Aqpred @ sN).reshape(3,1) # Sun pos body vector
    sBtrue = (Aqtrue @ sN).reshape(3,1)

    # Define sun sensor positions (one on each face)
    n1s = np.array([1, 0, 0]).reshape(-1,1)
    n2s = np.array([-1, 0, 0]).reshape(-1,1)
    n3s = np.array([0, 1, 0]).reshape(-1,1)
    n4s = np.array([0, -1, 0]).reshape(-1,1)
    n5s = np.array([0, 0, 1]).reshape(-1,1)
    n6s = np.array([0, 0, -1]).reshape(-1,1)

    Imax = 1

    p = np.dot(rvec.T, sN)
    if p > 0:
        illum = 1
    else:
        l = sqrt(np.linalg.norm(rvec)**2 - p**2)
        if l < 6371:
            illum = 0
        else:
            illum = 1
    
    if np.dot(n1s.T, sBtrue) < 0:
        Iss1m = 0
    else:
        ssnoise = rng.normal(0, 0.67*pi/180)
        Iss1m = illum*Imax*np.dot(n1s.T, sBtrue) + ssnoise
        Iss1p = illum*Imax*np.dot(n1s.T, sBest)
        yss1 = Iss1m - Iss1p
        Hss1 = np.hstack([n1s.T @ skew(sBest), np.zeros((1,3))])
        Rss1 = (0.67*pi/180)**2
        ylist.append(yss1)
        Hlist.append(Hss1)
        Rlist.append(Rss1)
    
    if np.dot(n2s.T, sBtrue) < 0:
        Iss2m = 0
    else:
        ssnoise = rng.normal(0, 0.67*pi/180)
        Iss2m = illum*Imax*np.dot(n2s.T, sBtrue) + ssnoise
        Iss2p = illum*Imax*np.dot(n2s.T, sBest)
        yss2 = Iss2m - Iss2p
        Hss2 = np.hstack([n2s.T @ skew(sBest), np.zeros((1,3))])
        Rss2 = (0.67*pi/180)**2
        ylist.append(yss2)
        Hlist.append(Hss2)
        Rlist.append(Rss2)

    if np.dot(n3s.T, sBtrue) < 0:
        Iss3m = 0
    else:
        ssnoise = rng.normal(0, 0.67*pi/180)
        Iss3m = illum*Imax*np.dot(n3s.T, sBtrue) + ssnoise
        Iss3p = illum*Imax*np.dot(n3s.T, sBest)
        yss3 = Iss3m - Iss3p
        Hss3 = np.hstack([n3s.T @ skew(sBest), np.zeros((1,3))])
        Rss3 = (0.67*pi/180)**2
        ylist.append(yss3)
        Hlist.append(Hss3)
        Rlist.append(Rss3)

    if np.dot(n4s.T, sBtrue) < 0:
        Iss4m = 0
    else:
        ssnoise = rng.normal(0, 0.67*pi/180)
        Iss4m = illum*Imax*np.dot(n4s.T, sBtrue) + ssnoise
        Iss4p = illum*Imax*np.dot(n4s.T, sBest)
        yss4 = Iss4m - Iss4p
        Hss4 = np.hstack([n4s.T @ skew(sBest), np.zeros((1,3))])
        Rss4 = (0.67*pi/180)**2
        ylist.append(yss4)
        Hlist.append(Hss4)
        Rlist.append(Rss4)

    if np.dot(n5s.T, sBtrue) < 0:
        Iss5m = 0
    else:
        ssnoise = rng.normal(0, 0.67*pi/180)
        Iss5m = illum*Imax*np.dot(n5s.T, sBtrue) + ssnoise
        Iss5p = illum*Imax*np.dot(n5s.T, sBest)
        yss5 = Iss5m - Iss5p
        Hss5 = np.hstack([n5s.T @ skew(sBest), np.zeros((1,3))])
        Rss5 = (0.67*pi/180)**2
        ylist.append(yss5)
        Hlist.append(Hss5)
        Rlist.append(Rss5)

    if np.dot(n6s.T, sBtrue) < 0:
        Iss6m = 0
    else:
        ssnoise = rng.normal(0, 0.67*pi/180)
        Iss6m = illum*Imax*np.dot(n6s.T, sBtrue) + ssnoise
        Iss6p = illum*Imax*np.dot(n6s.T, sBest)
        yss6 = Iss6m - Iss6p
        Hss6 = np.hstack([n6s.T @ skew(sBest), np.zeros((1,3))])
        Rss6 = (0.67*pi/180)**2
        ylist.append(yss6)
        Hlist.append(Hss6)
        Rlist.append(Rss6)

    # Magnetometer

    BECI = magfield(sol.t[j], rvec).reshape(3,1)
    BN = (BECI/np.linalg.norm(BECI)).reshape(3,1)
    BBtrue = (Aqtrue @ BN).reshape(3,1)
    mfnoise = 1.5 * pi/180

    BBmeas = BBtrue + rng.normal(0, mfnoise, (3,1))
    BBmeas /= np.linalg.norm(BBmeas)

    BBest = (Aqpred @ BN).reshape(3,1)
    ymf = BBmeas.reshape(3,1) - BBest.reshape(3,1)
    Hmf = np.hstack([skew(BBest), np.zeros((3, 3))])
    Rmf = (mfnoise**2) * np.eye(3)

    Hlist.append(Hmf)
    Rlist.append(Rmf)
    ylist.append(ymf)

    H = np.vstack(Hlist)
    R = scipy.linalg.block_diag(*Rlist)
    y = np.vstack(ylist)

    Pprev = P[:, :, j-1].reshape(6,6)
    wx, wy, wz = west.flatten()
    westsk = np.array([
    [ 0, -wz,  wy],
    [ wz,  0, -wx],
    [-wy,  wx,  0]
    ])
    F = np.zeros((6, 6))
    F[0:3, 0:3] = -westsk
    F[0:3, 3:6] = -np.eye(3)
    Phi = np.eye(6) + F * dt
    Pest = Phi @ Pprev @ Phi.T + Q_discrete

    S = H @ Pest @ H.T + R
    K = Pest @ H.T @ np.linalg.inv(S)

    dx = (K @ y).flatten()
    dth = dx[0:3]
    dbeta = dx[3:6]

    dqerror = np.array([1, 0.5*dth[0], 0.5*dth[1], 0.5*dth[2]]).reshape(-1,1)
    dqerror /= np.linalg.norm(dqerror)

    # Multiplicative Update
    qup = quatmul(xest[0:4], dqerror)
    #qup = quatmul(dqerror, xest)
    qup /= np.linalg.norm(qup)

    # Additive Update for gyro
    betaup = xest[4:7].flatten() + dbeta
    betaup = betaup.reshape(-1,1)

    I = np.eye(6)
    IKH = I - K @ H
    P[:, :, j] = IKH @ Pest @ IKH.T + K @ R @ K.T
    xhat[j, :] = np.concatenate([qup, betaup]).flatten()

    # Control Law

    westup = wmeas - betaup
    quateKF, weKF, RBLKF = errorfunc(sol.t[j], rvec, vvec, qup, westup, kep0)
    if quateKF[0] < 0:
        quateKF = -quateKF

    tauctrl = ctrllaw(quateKF, weKF)

    quateKFplt.append(quateKF.flatten())
    weKFplt.append(weKF.flatten())

quatetrueplt = np.array(quatetrueplt)
wetrueplt = np.array(wetrueplt)
hrwtrueplt = np.array(hrwtrueplt)
quateKFplt = np.array(quateKFplt)
weKFplt = np.array(weKFplt)
betahist = np.array(betahist)
#yplt = np.array(yplt)

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(time, quatetrueplt[:, 0], 'k--', label='Dotted = True Error')
plt.plot(time, quatetrueplt[:, 1], 'r--')
plt.plot(time, quatetrueplt[:, 2], 'g--')
plt.plot(time, quatetrueplt[:, 3], 'b--')
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(time[1:], quateKFplt[:, 1], 'r')
plt.plot(time[1:], quateKFplt[:, 2], 'g')
plt.plot(time[1:], quateKFplt[:, 3], 'b')
plt.plot(time[1:], quateKFplt[:, 0], 'k', label='Straight = MEKF error')
plt.xlabel("Time (seconds)")
plt.ylabel("Quaternion Error")
plt.legend()
plt.show(block=False)

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(time, wetrueplt[:, 0], 'k--', label='True Angular Velocity')
plt.plot(time, wetrueplt[:, 1], 'g--')
plt.plot(time, wetrueplt[:, 2], 'b--')
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(time[1:], weKFplt[:, 0], 'k', label = 'MEKF Angular Velocity Error')
plt.plot(time[1:], weKFplt[:, 1], 'g')
plt.plot(time[1:], weKFplt[:, 2], 'b')
plt.xlabel("Time (seconds)")
plt.ylabel("Angular Velocity Error")
plt.legend()
plt.show(block=False)

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(time[1:], betahist[:, 0], 'k', label = 'True bias')
plt.plot(time[1:], betahist[:, 1], 'b')
plt.plot(time[1:], betahist[:, 2], 'g')
plt.legend()
plt.subplot(2, 1, 2)
plt.plot(time, xhat[:, 4], 'k', label='Estimated Bias')
plt.plot(time, xhat[:, 5], 'b')
plt.plot(time, xhat[:, 6], 'g')
plt.xlabel("Time (seconds)")
plt.ylabel("Bias")
plt.legend()

plt.show()

#plt.figure()
#plt.plot(time[1:], yplt[:, 0])
#plt.plot(time[1:], yplt[:, 1])
#plt.plot(time[1:], yplt[:, 2])

#plt.show()