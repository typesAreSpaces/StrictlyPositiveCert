with(SolveTools, SemiAlgebraic, Parametric);

example1 := proc()
local i;
local f1 := 13465180504391033/20743227500000000*x + 46783161376177/163332500000000 + 47/8297291*x^6 - 244259/829729100*x^5 + 487010193/82972910000*x^4 - 460609659017/8297291000000*x^3 + 25233657328627/103716137500000*x^2;
local g1 := -2*x^6 + 5197/50*x^5 - 10361919/5000*x^4 + 9800205511/500000*x^3 - 536886326141/6250000*x^2 + 154852063736361/1250000000*x + 107498916600543/1250000000;

local sol := solve(
    {c > 0, c * f1 - g1 >= 0},
    {x}, 'parametric', 'real', 'parameters' = {c});

local _args := op(sol);
local last_condition := NULL;
local pos_coeff := 1;

    for i from 1 to nops(sol)/2 do
        if(evalb(_args[2*i] = [[x = x]])) then
            pos_coeff :=
            select(type,
                   [op(last_condition), op(evalf(_args[2*i - 1]))],
                   'numeric');
            pos_coeff := convert((pos_coeff[1] + pos_coeff[2])/2, rational, exact);
            break;
        end if;
        last_condition := evalf(simplify(_args[2*i - 1]));
    end do;

    print(pos_coeff);
    print(SemiAlgebraic([pos_coeff*f1 - g1 < 0], [x]));
end proc;

example2 := proc()
local i;
local f := x^2 + x + 1;

local sol := solve(
    {c > 0, c * f >= 0},
    {x}, 'parametric', 'real', 'parameters' = {c});

    print(sol);

local _args := op(sol);
local last_condition := NULL;
local pos_coeff := 1;

    for i from 1 to nops(sol)/2 do
        if(evalb(_args[2*i] = [[x = x]])) then
            pos_coeff :=
            select(type,
                   [op(last_condition), op(evalf(_args[2*i - 1]))],
                   'numeric');
            pos_coeff := convert((pos_coeff[1] + pos_coeff[2])/2, rational, exact);
            break;
        end if;
        last_condition := evalf(simplify(_args[2*i - 1]));
    end do;

    print(pos_coeff);
    print(SemiAlgebraic([pos_coeff*f < 0], [x]));
end proc;

#example1();
example2();
