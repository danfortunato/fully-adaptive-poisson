classdef Timer < handle

    properties
        depth = 0
        next_id = 0
        active = true
        event = emptyEvent()
        spool = emptyEvent()
    end

    methods ( Static )

        function obj = Timer()
            persistent localObj;
            if ( isempty(localObj) || ~isvalid(localObj) )
                localObj = obj;
            end
            obj = localObj;
        end

        function tic()
            timer = Timer();
            e = struct('tic', [], 'toc', [], 'message', [], 'depth', [], 'id', []);
            e.tic   = tic();
            e.depth = timer.depth;
            e.id    = timer.next_id;
            timer.depth   = timer.depth   + 1;
            timer.next_id = timer.next_id + 1;
            timer.event(end+1) = e;
        end

        function toc(message)
            timer = Timer();

            if ( isempty(timer.event) )
                error('toc() has no matching tic().');
            end

            e = timer.event(end);
            e.toc = toc(e.tic);
            e.message = message;
            timer.spool(end+1) = e;
            timer.event(end) = [];
            timer.depth = timer.depth - 1;
            if ( timer.active && timer.depth == 0 )
                timer.next_id = 0;
                Timer.print();
            end
        end

        function log(message)
            timer = Timer();
            e = struct('tic', [], 'toc', [], 'message', [], 'depth', [], 'id', []);
            e.message = message;
            e.depth = timer.depth;
            e.id    = timer.next_id;
            timer.next_id = timer.next_id + 1;
            timer.spool(end+1) = e;
            if ( timer.active && timer.depth == 0 )
                timer.next_id = 0;
                Timer.print();
            end
        end

        function reset()
            timer = Timer();
            timer.depth = 0;
            timer.next_id = 0;
            timer.active = true;
            timer.event = emptyEvent();
            timer.spool = emptyEvent();
        end

        function on()
            timer = Timer();
            timer.active = true;
        end

        function off()
            timer = Timer();
            timer.active = false;
        end

    end

    methods ( Static, Access = private )

        function print()
            timer = Timer();
            if ( timer.active )
                [~, idx] = sort([timer.spool.id]);
                for k = 1:length(timer.spool)
                    e = timer.spool(idx(k));
                    pre = '';
                    if ( e.depth > 0 )
                        pre = ['  ' repmat('    ', 1, e.depth-1) '─ '];
                    end
                    if ( isempty(e.toc) )
                        fprintf([pre e.message '\n']);
                    else
                        fprintf([pre e.message ': %gs\n'], e.toc);
                    end
                end
                timer.spool(1:end) = [];
            end
        end

    end

end

function event = emptyEvent
    persistent emptyEvent
    if ( isempty(emptyEvent) )
        fields = {'tic', 'toc', 'message', 'depth', 'id'};
        fields{2,1} = {};
        emptyEvent = struct(fields{:});
    end
    event = emptyEvent;
end
