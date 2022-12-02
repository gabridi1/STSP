using Gurobi, JuMP, Graphs, GraphPlot, Random

#-------------------------Funções-----------------------------------

cortes = 0

function menorCiclo(x, n)
    shortest_subtour = collect(1: n)
    unvisited = Set(collect(1: n))

    while !isempty(unvisited)
        vizinhos = unvisited
        this_cycle = Int[]
        while !isempty(vizinhos)
            current = pop!(vizinhos)
            push!(this_cycle, current)

            if length(this_cycle) > 1
                delete!(unvisited, current)
            end

            vizinhos = [j for j in 1: n if j in unvisited && value(x[current, j]) > 2/3]
        end

        if length(this_cycle) < length(shortest_subtour)
            shortest_subtour = this_cycle
        end
    end
    
    return shortest_subtour
end

function generate_distance_matrix(n; random_seed = 1)
    rng = Random.MersenneTwister(random_seed)
    X = 100 * rand(rng, n)
    Y = 100 * rand(rng, n)
    d = [sqrt((X[i] - X[j])^2 + (Y[i] - Y[j])^2) for i in 1: n, j in 1: n]
    return X, Y, d
end

function plotGraph(graph)
    nodelabel = 1:nv(graph)
    p = gplot(graph, nodelabel = nodelabel)
    display(p)
end

function my_callback_function_lazy(cb_data)
    x_val = callback_value.(cb_data, x)

    status = MOI.get(model, MOI.CallbackNodeStatus(cb_data))

    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        shortest_subtour = menorCiclo(x_val, n)
        if 1 < length(shortest_subtour) < n
            con = @build_constraint(sum(x[i, j] for i in 1: n, j in 1: n if i in shortest_subtour && j in shortest_subtour && j < i) <= length(shortest_subtour) - 1)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con) 
            global cortes += 1    
            return 
        end
    end
end

#---------------------------------Main--------------------------------------

n = 100

X, Y, d = generate_distance_matrix(n)

model = Model(Gurobi.Optimizer)

set_optimizer_attribute(model, "LazyConstraints", 1)
set_optimizer_attribute(model, "Cuts", 0)
set_optimizer_attribute(model, "Presolve", 0)
set_optimizer_attribute(model, "Heuristics", 0)
set_optimizer_attribute(model, "OutputFlag", 0)

@variable(model, x[1: n, 1: n], Bin, Symmetric)
@objective(model, Min, sum(d .* x) / 2)
@constraint(model, [i in 1: n], sum(x[i, :]) == 2)
@constraint(model, [i in 1: n], x[i, i] == 0)

MOI.set(model, MOI.LazyConstraintCallback(), my_callback_function_lazy)
optimize!(model)
time = solve_time(model)
nos = node_count(model)

G = SimpleGraph(n)

for i in 1: n
    for j in 1: n
        if j < i
            if value(x[i, j]) > 0.5
                add_edge!(G, i, j)
                println("Arco ($i, $j)")
            end
        end 
    end 
end

plotGraph(G)
println("z = ", objective_value(model))
println("Nós = ", nos)
println("Cortes = ", cortes)
print("Tempo = ", time)
