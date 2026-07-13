classdef TimeStepper
    %TIMESTEPPER Deterministic fixed-step time-loop helper.
    methods (Static)
        function state = ForwardEuler(rhs, state, dt, steps)
            for k = 1:steps
                state = state + dt * rhs(state, k * dt);
            end
        end
    end
end
