using ModelingToolkit, ModelingToolkitStandardLibrary

@parameters t
const Rot = ModelingToolkitStandardLibrary.Mechanical.Rotational
const B = ModelingToolkitStandardLibrary.Blocks

@component function LinearSpacecraftAttitude(
    ; name, Jx=100.0, Jy=100.0, Jz=100.0, u0=[0.0,0.0,0.0], ω0=[0.0,0.0,0.0], ω̇0=[0.0,0.0,0.0]
)

    @named Ix = Rot.Inertia(J=Jx, phi_start=u0[1], w_start=ω0[1], a_start=ω̇0[1])
    @named Iy = Rot.Inertia(J=Jy, phi_start=u0[2], w_start=ω0[2], a_start=ω̇0[2])
    @named Iz = Rot.Inertia(J=Jz, phi_start=u0[3], w_start=ω0[3], a_start=ω̇0[3])

    @named x_flange_a = Rot.Flange()
    @named y_flange_a = Rot.Flange()
    @named z_flange_a = Rot.Flange()

    @named x_flange_b = Rot.Flange()
    @named y_flange_b = Rot.Flange()
    @named z_flange_b = Rot.Flange()

    @named ϕ_sensor = Rot.AngleSensor()
    @named θ_sensor = Rot.AngleSensor()
    @named ψ_sensor = Rot.AngleSensor()

    ps = @parameters Jx=Jx Jy=Jy Jz=Jz u0=u0 ω0=ω0 ω̇0=ω̇0

    D = Differential(t)

    eqs = [
        connect(x_flange_a, Ix.flange_a),
        connect(y_flange_a, Iy.flange_a),
        connect(z_flange_a, Iz.flange_a),

        connect(Ix.flange_b, x_flange_b),
        connect(Iy.flange_b, y_flange_b),
        connect(Iz.flange_b, z_flange_b),

        connect(x_flange_b, ϕ_sensor.flange),
        connect(y_flange_b, θ_sensor.flange),
        connect(z_flange_b, ψ_sensor.flange),
    ]

    compose(
        ODESystem(eqs, t, [], ps; name = name),
        Ix, Iy, Iz,
        x_flange_a, y_flange_a, z_flange_a,
        x_flange_b, y_flange_b, z_flange_b,
        ϕ_sensor, θ_sensor, ψ_sensor
    )
end

@component function SpacecraftAttitude(
    ; name, Jx=100.0, Jy=100.0, Jz=100.0, u0=zeros(3), ω0=zeros(3), ω̇0=zeros(3)
)

    @named Mx = B.RealInput()
    @named My = B.RealInput()
    @named Mz = B.RealInput()

    @named phi_x = B.RealOutput()
    @named phi_y = B.RealOutput()
    @named phi_z = B.RealOutput()

    sts = @variables ϕ(t)=u0[1] θ(t)=u0[2] ψ(t)=u0[3] ωx(t)=ω0[1] ωy(t)=ω0[2] ωz(t)=ω0[3] ω̇x(t)=ω̇0[1] ω̇y(t)=ω̇0[2] ω̇z(t)=ω̇0[3]

    ps = @parameters Jx=Jx Jy=Jy Jz=Jz u0=u0 ω0=ω0 ω̇0=ω̇0

    D = Differential(t)

    eqs = [
        phi_x.u ~ ϕ,
        phi_y.u ~ θ,
        phi_z.u ~ ψ,

        D(ϕ) ~ ωx + ωz * tan(θ)*cos(ϕ) + ωy*tan(θ)*sin(ϕ),
        D(θ) ~ ωy*cos(ϕ) - ωz*sin(ϕ),
        D(ψ) ~ ωz*sec(θ)*cos(ϕ) + ωy*sec(θ)*sin(ϕ),

        D(ωx) ~ ω̇x,
        D(ωy) ~ ω̇y,
        D(ωz) ~ ω̇z,

        Jx * ω̇x ~ Mx.u + (Jy - Jz)*ωy*ωz,
        Jy * ω̇y ~ My.u + (Jz - Jx)*ωx*ωz,
        Jz * ω̇z ~ Mz.u + (Jx - Jy)*ωx*ωy,
    ]

    compose(
        ODESystem(eqs, t, sts, ps; name = name), Mx, My, Mz, phi_x, phi_y, phi_z
    )
end