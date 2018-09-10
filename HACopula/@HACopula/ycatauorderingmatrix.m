function tauOrdMat = ycatauorderingmatrix(obj)
%YCATAUORDERINGMATRIX - TauOrdering of the yca for each pair of leaves
%
%
% Copyright 2018 Jan Gorecki

tauOrdMat = zeros(obj.Dim, obj.Dim);
tauOrdMat = ycatauorderingmatrixrec(obj, tauOrdMat);
tauOrdMat = tauOrdMat + tauOrdMat'; % make it symmetric

end

function tauOrdMat = ycatauorderingmatrixrec(obj, tauOrdMat)
    for i = 1:length(obj.Child)
        for j = (i+1):length(obj.Child)
            if isa(obj.Child{i}, 'HACopula')
                iChildren = obj.Child{i}.Leaves;
            else
                iChildren = obj.Child{i};
            end
            if isa(obj.Child{j}, 'HACopula')
                jChildren = obj.Child{j}.Leaves;
            else
                jChildren = obj.Child{j};
            end
            tauOrdMat(iChildren, jChildren) = obj.TauOrdering;
        end
    end
    for i = 1:length(obj.Child)
        if isa(obj.Child{i}, 'HACopula')
            tauOrdMat = ycatauorderingmatrixrec(obj.Child{i}, tauOrdMat);
        end
    end
end
