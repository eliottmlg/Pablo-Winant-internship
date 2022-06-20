
struct Iterationlog
    header::Array{String, 1}
    keywords::Array{Symbol, 1}
    types::Vector
    entries::Array{Any, 1}
    fmt::String
end

# function Iterationlog()
#     header = [ "It", "ϵₙ", "ηₙ=|xₙ-xₙ₋₁|", "λₙ=ηₙ/ηₙ₋₁", "Time", "Newton steps"]
#     keywords = [:it, :err, :gain, :time, :nit]
#     Iterationlog(header, keywords, [])
# end

function Iterationlog(;kw...)
    keywords = Symbol[]
    header = String[]
    elts = []
    types= []
    for (k,(v,t)) in kw
        push!(types, t)
        push!(keywords, k)
        push!(header, v)
        if t<:Int
            push!(elts, "{:>8}")
        elseif t<:Real
            push!(elts, "{:>12.4e}")
        else
            error("Don't know how to format type $t")
        end
    end
    fmt = join(elts, " | ")
    Iterationlog(header, keywords, types, [], fmt)
end

function append!(log::Iterationlog; verbose=true, entry...)
    d = Dict(entry)
    push!(log.entries, d)
    verbose && show_entry(log, d)
end

function initialize(log::Iterationlog; verbose=true, message=nothing)
    verbose && show_start(log; message=message)
end

function finalize(log::Iterationlog; verbose=true)
    verbose && show_end(log)
end

function show(log::Iterationlog)
    show_start(log)
    for entry in log.entries
        show_entry(log,entry)
    end
    show_end(log)
end

function show_start(log::Iterationlog; message=nothing)
    # println(repeat("-", 66))
    # @printf "%-6s%-16s%-16s%-16s%-16s%-5s\n" "It" "ϵₙ" "ηₙ=|xₙ-xₙ₋₁|" "λₙ=ηₙ/ηₙ₋₁" "Time" "Newton steps"
    # println(repeat("-", 66))

    l = []
    for t in log.types
        if t<:Int
            push!(l, "{:>8}")
        elseif t<:Real
            push!(l, "{:>12}")
        end
    end

    a = join(l, " | ")
    s = format(a, log.header...)

    if message !== nothing
        println(repeat("-", length(s)))
        println(message)
    end
    println(repeat("-", length(s)))

    println(s)
    println(repeat("-", length(s)))
end


function show_entry(log::Iterationlog, entry)

    vals = [entry[k] for k in log.keywords]
    printfmt(log.fmt, vals...)
    print("\n")

end

function show_end(log::Iterationlog)
    println(repeat("-", 66))
end

###

mutable struct IterationTrace
    trace::Array{Any,1}
end


####################### Same beginning as Time iteration ##############################

