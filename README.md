# TOV_MI_Tidal
# Anisotropic CGC Dark Stars in Rastall Gravity

This repository contains programs that analyze the following compact star model:
- **EoS**: Simplified general Chaplygin dark star (input parameters A & B)
- **Anisotropy factor**: model by Horvat et al. model (input parameter alpha)
- **Modified gravity framework**: Rastall gravity (input parameter beta)

The programs use the RK4 method to compute the following properties:
- TOV equations,
- Moment of inertia
- Electric tidal deformability

Programs are available in Fortran (.f) and Python notebook (.ipynb)

## Fortran program: TOVMITIDAL_Rastall_anisotropic.f
Modify anisotropic and Rastall parameters and PCC limits here:

```fortran
c---  !change anisotropy, rastall, PCC parameters here: 
      alp = 0.0                 ! anisotropy factor (isotropic = 0)
      rast = 0.05               ! Rastall parameter (GR = 0)
      PCCmax = 500              ! max PCC to stop printing
      PCCprof = 101.D0          ! PCC point to save star profile data 
```

Modify Chaplygin EoS and its parameters inside the `FED` function:
```fortran
c---  !change Chaplygin dark star parameters here
      A = 0.3D0                     ! dark star CG parameter A=0.3
      B = 6.0D-20 * (conv**2)       ! dark star CG parameter B=6e-20
```

The program automatically stops once star structure disobeys energy conditions.

## Python notebook: TOV_MI_Tidal_Rastall.ipynb
Modify anisotropic and Rastall parameters here:
```python
ani = 0     # set 0 for isotropes
ras = 0     # set 0 for GR
```
Modify Chaplygin EoS and its parameters inside the `eden` function:
```python
    def eden(self, P):
        conv = 7.5534e11
        a = 0.3
        b = 6.0 * 1e-20 * conv * conv
```
