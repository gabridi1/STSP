using Graphs, GraphPlot, Multigraphs, GLPK, JuMP, Random


#----------------------Constantes-------------------------

const K = 1000    #Iterações
const c1 = 2      #Descida exponencial
const c2 = 50     #Divisor
const c3 = 0.00001  #Teta
const c4 = 0.001  #Gap

@time begin
    
#------------------------Funções----------------------------

function plotGraph(graph)
    nodelabel = 1: nv(graph)
    p = gplot(graph, nodelabel = nodelabel)
    display(p)
end

function generate_distance_matrix(n; random_seed = 1)
    rng = Random.MersenneTwister(random_seed)
    X = 100 * rand(rng, n)
    Y = 100 * rand(rng, n)
    d = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1: n, j in 1: n]
    return X, Y, d
end

function custosLagrangeanos(c, u, n)
    cl = zeros(n, n)

    for i in 1: n, j in 1: n
        cl[i, j] = c[i, j] - u[i] - u[j]
    end 
    
    return cl 
end

function umaArvore(cl, n)
    G = CompleteGraph(n - 1)
    T_out = prim_mst(G, cl[2: n, 2: n])
    T = SimpleGraph(nv(G))

    for edge in T_out
        add_edge!(T, edge)
    end

    um_vizinhos = sortperm(cl[1, 2: n])

    uma_arvore = SimpleGraph(n)

    add_edge!(uma_arvore, 1, 1 + um_vizinhos[1])
    add_edge!(uma_arvore, 1, 1 + um_vizinhos[2])
    
    for i in 1: n - 1, j in 1: n - 1
        if has_edge(T, i, j)
            add_edge!(uma_arvore, i + 1, j + 1)
        end
    end
       
    return uma_arvore  
end   

function emparelhamento(cm, n)
    model = Model(GLPK.Optimizer)

    @variable(model, x[1: n, 1: n], Bin, Symmetric)
    @constraint(model, linha[i in 1: n], sum(x[i, j] for j in 1: n) == 1)
    @constraint(model, diagonal[i in 1: n], x[i, i] == 0)
    @objective(model, Min, sum(cm[i, j] * x[i, j] for i in 1: n, j in 1: n))

    optimize!(model)

    return x
end

function fatorDois(c_comp, n)
    #1 - Árvore ótima e duplicar
    G = CompleteGraph(n)
    T_out = prim_mst(G, c_comp)
    T = Multigraph(nv(G))

    for edge in T_out
        add_edge!(T, edge)
        add_edge!(T, edge)
    end

    #2 -Euleriano
    euleriano = Array{Int64}(undef, 2 * n - 1)
    euleriano[1] = 1
    current1 = 1
    current2 = 1

    
    for i in 1: 2 * n - 2
        current1 = current2
        vizinhos = [j for j in 1: n if j in 1: n && has_edge(T, current1, j)]
        
        controle = true
        while controle
            if length(vizinhos) > 0
                current2 = pop!(vizinhos)
            end

            rem_edge!(T, current1, current2)

            v = connected_components(T)
            
            cc = 0
            for i in 1: length(v)
                if length(v[i]) > 1
                    cc += 1
                end
            end

            if cc > 1
                add_edge!(T, current1, current2)
            else
                controle = false
            end
        end
        
        euleriano[i + 1] = current2
    end
    
    #3 - Hamiltoniano
    euleriano = unique(euleriano)
    H = SimpleGraph(n)
    
    for i in 1: length(euleriano) - 1
        add_edge!(H, euleriano[i], euleriano[i + 1])
    end
    
    add_edge!(H, 1, euleriano[length(euleriano)])

    return H   
end

function christofides(c_comp, n)
    
    #1 - Árvore ótima
    G = CompleteGraph(n)
    T_out = prim_mst(G, c_comp)
    T = Multigraph(nv(G))

    for edge in T_out
        add_edge!(T, edge)
    end

    #2 - Emparelhamento
    v = Array{Int64}(undef, n)
    cont = 1
    for i in 1: n
        if length(neighbors(T, i)) % 2 != 0
            v[cont] = i
            cont += 1
        end 
    end
    
    v = v[1: cont - 1]
    
    if(cont > 1)
        cm = c_comp[v, v]  
        
        x = emparelhamento(cm, length(v))
        
        na = n - 1
        for i in 1: length(v), j in 1: length(v)
            if value(x[i, j]) > 0.5
                if j < i
                    add_edge!(T, v[i], v[j])
                    na += 1
                end
            end
        end
    end

    #3 -Euleriano
    euleriano = Array{Int64}(undef, na + 1)
    euleriano[1] = 1
    current1 = 1
    current2 = 1

    for i in 1: na
        current1 = current2
        vizinhos = [j for j in 1: n if j in 1:n && has_edge(T, current1, j)]
        
        controle = true
        while controle
            if length(vizinhos) > 0
                current2 = pop!(vizinhos)
            end

            rem_edge!(T, current1, current2)

            v = connected_components(T)       
            cc = 0

            for i in 1: length(v)
                if length(v[i]) > 1
                    cc += 1
                end
            end

            if cc > 1
                add_edge!(T, current1, current2)
            else
                controle = false
            end
        end
        
        euleriano[i + 1] = current2
    end
    
    #4 - Hamiltoniano
    hamiltoniano = unique(euleriano)
    H = SimpleGraph(n)
    
    for i in 1: length(hamiltoniano) - 1
        add_edge!(H, hamiltoniano[i], hamiltoniano[i + 1])
    end
    
    add_edge!(H, 1, hamiltoniano[length(hamiltoniano)])
  
    return H   
end

function limiteDual(c, u, n)
    cl = custosLagrangeanos(c, u, n)
    uma_arvore = umaArvore(cl, n)
    
    zd = 2 * sum(u[i] for i in 1: n)
    
    for i in 1: n, j in 1: n
        if j < i
            if has_edge(uma_arvore, i, j) 
                zd += cl[i, j]
            end
        end
    end
        
    return zd, uma_arvore
end

function limitePrimal(c, uma_arvore, n, u)
    c_comp = Array{Float64, 2}(undef, n, n)
    
    for i in 1: n, j in 1: n
        if has_edge(uma_arvore, i, j) 
            c_comp[i, j] = 0
        else
            c_comp[i, j] = c[i, j]
        end
    end

    cl = custosLagrangeanos(c, u, n)

    H = christofides(c_comp, n)
    #H = fatorDois(c_comp, n)
    #H = christofides(cl, n)
    #H = fatorDois(cl, n)
    #H = christofides(c, n)
    #H = fatorDois(c, n)
    
    zp = 0 
    for i in 1: n, j in 1: n
        if j < i
            if has_edge(H, i, j) 
                zp += c[i, j]
            end
        end
    end

    return zp
end

function relaxacaoLagrangeana(c, n)
    max_iter = K

    z_lp = Inf
    z_ld = - Inf
    
    u = zeros(n)
    
    contador = 0 
    theta = 2
        
    #z_lp = limitePrimal(c, 1, n, u)
    
    for k = 1: max_iter
        zd, uma_arvore = limiteDual(c, u, n)
        zp = limitePrimal(c, uma_arvore, n, u)
        
        if zp < z_lp
            z_lp = zp
        end

        if zd > z_ld
            z_ld = zd
            contador = 0
        else
            contador += 1
        end
        
        #println(z_lp, "   ", z_ld)
        
        if contador == max_iter/c2
            theta /= c1
            contador = 0
        end
            
        if theta <= c3
            println("Saindo sem otimalidade!")
            println("Critério: theta pequeno")
            break
        end
        
        dif = zeros(n)
        for i in 1: n
            dif[i] = 2 - length(neighbors(uma_arvore, i))
        end

        soma = sum(dif[i]^2 for i in 1: n)
        
        if soma != 0
            mu = theta * (z_lp - zd) / soma
        else
            mu = 0
            println("Saindo por otimalidade!")
            println("Critério: Dual otimalidade")
            z_lp = zd
            break
        end        
    
        u += mu * dif

        opt_gap = (z_lp - z_ld) / z_lp
        if opt_gap < c4
            println("Saindo por otimalidade!")
            println("Critério: gap de otimalidade")
            break
        end
        
        if k == max_iter
            println("Saindo sem otimalidade!")
            println("Critério: número de iterações")
        end
    end
          
    return z_lp, z_ld
end
    
#-----------------------------Main----------------------------

n = 100
    
X, Y, d = generate_distance_matrix(n)
    
z_lp, z_ld = relaxacaoLagrangeana(d, n)
    
end