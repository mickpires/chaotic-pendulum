using Plots

mutable struct pendulo
    l::Number
    h::Number
    q::Number # dragging coeficient
    Ω::Number #dragging frequency
    F::Number
    periodo::Number
    tempo_final::Number
    num_pontos::Int64
    θ::Vector{Number}
    ω::Vector{Number}
    Δt::Vector{Number}

    function pendulo(;initial_θ = 0,h::Number)
        l = 9.81
        q = 1/3
        Ω = 2/3
        F = 1.4
        g = 9.81
        periodo = 2 * π / sqrt(g/l)
        tempo_final = 10 * periodo
        num_pontos = ceil(Int64,tempo_final/h)
        θ = zeros(num_pontos)
        θ[1] = initial_θ
        ω = zeros(num_pontos)
        Δt = range(0,tempo_final,num_pontos)
        new(l,h,q,Ω,F,periodo,tempo_final,num_pontos,θ,ω,Δt)
    end
end

function euler(p::pendulo)
    g = 9.81
    for i = 1:p.num_pontos-1
        p.ω[i+1] = p.ω[i] + (-g * sin(p.θ[i])/ p.l - p.q * p.ω[i] + p.F * sin(p.Ω * p.Δt[i])) * p.h
        p.θ[i+1] = p.θ[i] + p.ω[i+1] * p.h
    end
    return nothing
end
h = 1e-3
pendulo1 = pendulo(initial_θ = 0.2,h=h)
euler(pendulo1)
pendulo2 = pendulo(initial_θ=0.21,h=h)
euler(pendulo2)

plt=plot(pendulo1.Δt, [pendulo1.θ,pendulo2.θ], 
label=["pendulo 1" "pendulo 2"], 
linewidth=2,
xlabel="Tempo (s)",
ylabel="θ")
display(plt)

espaço_fasico=plot([pendulo1.θ,pendulo2.θ],[pendulo1.ω,pendulo2.ω],
label=["Pêndulo 1" "Pêndulo 2"],
linewidth=2,
xlabel="θ",
ylabel="ω",
title="Espaço fásico dos dois pêndulos")