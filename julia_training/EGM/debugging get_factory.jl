
using Pkg
using Interpolations
using Plots
using Dolo
using Dolo: UNormal, discretize 
using StaticArrays
using QuantEcon
using LinearAlgebra
import Dolang

# wrong arguments, get_factory, gen_gufun, RECIPES #

## TO RUN BEFORE THE REST 
filename = "C:/Users/t480/GitHub/Pablo-Winant-internship/Dolo.jl/examples/models/consumption_savings_iid.yaml"
readlines(filename)
CSIID = yaml_import(filename)
_symbol_dict(x) = x
_symbol_dict(d::AbstractDict) =
    Dict{Symbol,Any}([(Symbol(k), _symbol_dict(v)) for (k, v) in d])
const src_path = dirname(@__FILE__)
const pkg_path = dirname(src_path)
Pkg.add("YAML")
import YAML; using YAML: load_file, load
const RECIPES = _symbol_dict(load_file(joinpath(src_path, "recipes.yaml")))
specs = RECIPES[:dtcc][:specs]
specs[Symbol("half_transition")]
specs[Symbol("half_transition")][Symbol("eqs")]
keys(specs)

# I run only the part that deals with argument in get_factory()
Symbol("half_transition") in keys(specs) 
Dolo.get_factory(CSIID, "half_transition") 
recipe = specs[Symbol("half_transition")]
symbols = Dolo.get_symbols(CSIID)
arguments = Dolo.OrderedDict(
            Symbol(l[3]) => [Dolang.stringify(e,l[2]) for e in symbols[Symbol(l[1])]]
            for l in recipe[:eqs] if !(l[1]=="parameters")
        )
# No problem it returns the right arguments, so it must be above in the function
# I am calling the function from the wrong file 

# If I run the entire function get_factory(model::Model, eq_type::String) instead:
function get_factory(model::Model, eq_type::String)
    if eq_type == "arbitrage"
        defs_0 = get_definitions(model; stringify=true)
        defs_1 = get_definitions(model; tshift=1, stringify=true)
        definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in merge(defs_0, defs_1)])

        eqs, eq_lb, eq_ub = get_equation_block(model, eq_type)

        symbols = get_symbols(model)

        #TODO : fix crazy bug: it doesn't work without the trailing underscore !
        equations = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
            Symbol(string("out_", i, "_")) => eqs[i]
            for i =1:length(eqs)
        )
        arguments = OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :x => [stringify(e,0) for e in symbols[:controls]],
            :M => [stringify(e,1) for e in symbols[:exogenous]],
            :S => [stringify(e,1) for e in symbols[:states]],
            :X => [stringify(e,1) for e in symbols[:controls]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
        
        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))

        equations_lb = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i, "_")) => eq_lb[i]
                for i =1:length(eqs)
            )
        equations_ub = OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i,"_")) => eq_ub[i]
                for i =1:length(eqs)
            )
        arguments_bounds = OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )

        # definitions: we should remove definitions depending on controls
        # definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in defs_0] )
        definitions = OrderedDict{Symbol, SymExpr}()
        ff_lb = FunctionFactory(equations_lb, arguments_bounds, definitions, :controls_lb)
        ff_ub = FunctionFactory(equations_ub, arguments_bounds, definitions, :controls_ub)

        return (ff, ff_lb, ff_ub)
           
    # elseif eq_type=="transition"
    #     # defs_0  = get_definitions(model; stringify=true)
    #     defs_m1 = get_definitions(model; tshift=-1, stringify=true)

    #     definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

    #     equations = get_assignment_block(model, eq_type)
    #     symbols = get_symbols(model)
    #     arguments = OrderedDict(
    #         :m => [stringify(e,-1) for e in symbols[:exogenous]],
    #         :s => [stringify(e,-1) for e in symbols[:states]],
    #         :x => [stringify(e,-1) for e in symbols[:controls]],
    #         :M => [stringify(e,0) for e in symbols[:exogenous]],
    #         
    #     )
    #     ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
    #     return ff
    else
        specs = RECIPES[:dtcc][:specs]
        if !(Symbol(eq_type) in keys(specs))
            error("No spec for equation type '$eq_type'")
        end
        recipe = specs[Symbol(eq_type)]
        if !(:target in keys(recipe))
            error("Not implemented")
        end
        equations = get_assignment_block(model, eq_type)
        symbols = get_symbols(model)
        times = unique([e[2] for e in recipe[:eqs]])
        defs = [get_definitions(model; tshift=t, stringify=true) for t in times]
        definitions = OrderedDict{Symbol, SymExpr}(
            [ (Dolang.stringify(k), v) for (k,v) in merge(defs...)]
        )
        arguments = OrderedDict(
            Symbol(l[3]) => [stringify(e,l[2]) for e in symbols[Symbol(l[1])]]
            for l in recipe[:eqs] if !(l[1]=="parameters")
        )
        arguments[:p] = [stringify(e) for e in symbols[:parameters]]

        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
    end

end
# I get:
F_g = Dolo.get_factory(CSIID, "half_transition")
# moving specs = RECIPES[:dtcc][:specs] outside if-loop
F_g = Dolo.get_factory(CSIID, "half_transition")
# still... even though the following is true
(Symbol("half_transition") in keys(specs))
# since running the function alone generates the error, it is not necessarily 
# a problem of using the wrong file. It must be in the function itself.

# I rewrite get_factory() for each equation
function get_factory_noloop(model::Model, eq_type::String)
    if eq_type == "arbitrage"
        defs_0 = Dolo.get_definitions(model; stringify=true)
        defs_1 = Dolo.get_definitions(model; tshift=1, stringify=true)
        definitions = Dolo.OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in merge(defs_0, defs_1)])

        eqs, eq_lb, eq_ub = get_equation_block(model, eq_type)

        symbols = get_symbols(model)

        #TODO : fix crazy bug: it doesn't work without the trailing underscore !
        equations = Dolo.OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
            Symbol(string("out_", i, "_")) => eqs[i]
            for i =1:length(eqs)
        )
        arguments = Dolo.OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :x => [stringify(e,0) for e in symbols[:controls]],
            :M => [stringify(e,1) for e in symbols[:exogenous]],
            :S => [stringify(e,1) for e in symbols[:states]],
            :X => [stringify(e,1) for e in symbols[:controls]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
        
        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))

        equations_lb = Dolo.OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i, "_")) => eq_lb[i]
                for i =1:length(eqs)
            )
        equations_ub = Dolo.OrderedDict{Symbol, Union{Expr, Number, Symbol}}(
                Symbol(string("out_", i,"_")) => eq_ub[i]
                for i =1:length(eqs)
            )
        arguments_bounds = Dolo.OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )

        # definitions: we should remove definitions depending on controls
        # definitions = OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in defs_0] )
        definitions = Dolo.OrderedDict{Symbol, SymExpr}()
        ff_lb = FunctionFactory(equations_lb, arguments_bounds, definitions, :controls_lb)
        ff_ub = FunctionFactory(equations_ub, arguments_bounds, definitions, :controls_ub)

        return (ff, ff_lb, ff_ub)
           
    elseif eq_type=="transition"
         # defs_0  = get_definitions(model; stringify=true)
        defs_m1 = Dolo.get_definitions(model; tshift=-1, stringify=true)

        definitions = Dolo.OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

        equations = get_assignment_block(model, eq_type)
        symbols = get_symbols(model)
        arguments = Dolo.OrderedDict(
             :m => [stringify(e,-1) for e in symbols[:exogenous]],
             :s => [stringify(e,-1) for e in symbols[:states]],
             :x => [stringify(e,-1) for e in symbols[:controls]],
             :M => [stringify(e,0) for e in symbols[:exogenous]],
             :p => [stringify(e) for e in symbols[:parameters]]
             )
             
        ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
        return ff
    elseif eq_type=="expectation"
        # defs_0  = get_definitions(model; stringify=true)
       defs_m1 = Dolo.get_definitions(model; tshift=-1, stringify=true)

       definitions = Dolo.OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

       equations = get_assignment_block(model, eq_type)
       symbols = get_symbols(model)
       arguments = Dolo.OrderedDict(
            :M => [stringify(e,1) for e in symbols[:exogenous]],
            :S => [stringify(e,1) for e in symbols[:states]],
            :X => [stringify(e,1) for e in symbols[:controls]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
            
       ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
       return ff
    elseif eq_type=="half_transition"
        # defs_0  = get_definitions(model; stringify=true)
       defs_m1 = Dolo.get_definitions(model; tshift=-1, stringify=true)

       definitions = Dolo.OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

       equations = get_assignment_block(model, eq_type)
       symbols = get_symbols(model)
       arguments = Dolo.OrderedDict(
            :m => [stringify(e,-1) for e in symbols[:exogenous]],
            :a => [stringify(e,-1) for e in symbols[:poststates]],
            :M => [stringify(e,0) for e in symbols[:exogenous]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
            
       ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
       return ff
    elseif eq_type=="direct_response_egm"
        # defs_0  = get_definitions(model; stringify=true)
       defs_m1 = Dolo.get_definitions(model; tshift=-1, stringify=true)

       definitions = Dolo.OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

       equations = get_assignment_block(model, eq_type)
       symbols = get_symbols(model)
       arguments = Dolo.OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :a => [stringify(e,0) for e in symbols[:poststates]],
            :z => [stringify(e,0) for e in symbols[:expectations]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
            
       ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
       return ff
    elseif eq_type=="reverse_state"
        # defs_0  = get_definitions(model; stringify=true)
       defs_m1 = Dolo.get_definitions(model; tshift=-1, stringify=true)

       definitions = Dolo.OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

       equations = get_assignment_block(model, eq_type)
       symbols = get_symbols(model)
       arguments = Dolo.OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :a => [stringify(e,0) for e in symbols[:poststates]],
            :x => [stringify(e,0) for e in symbols[:controls]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
            
       ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
       return ff
    elseif eq_type=="direct_response"
        # defs_0  = get_definitions(model; stringify=true)
       defs_m1 = Dolo.get_definitions(model; tshift=-1, stringify=true)

       definitions = Dolo.OrderedDict{Symbol, SymExpr}([ (Dolang.stringify(k), v) for (k,v) in  defs_m1 ])

       equations = get_assignment_block(model, eq_type)
       symbols = get_symbols(model)
       arguments = Dolo.OrderedDict(
            :m => [stringify(e,0) for e in symbols[:exogenous]],
            :s => [stringify(e,0) for e in symbols[:states]],
            :z => [stringify(e,0) for e in symbols[:expectations]],
            :p => [stringify(e) for e in symbols[:parameters]]
            )
            
       ff = FunctionFactory(equations, arguments, definitions, Symbol(eq_type))
       return ff
    end

end
# I get:
F_g = get_factory_noloop(CSIID, "half_transition")
# so I replace get_definitions by Dolo.get_definitions and rerun 

code_g = Dolang.gen_generated_gufun(F_g)
