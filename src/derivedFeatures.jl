function paddedDiffsqr!(x, y)
    #return [0 diff(x,2)].^2
    for i=1:length(y)-1
        d = (x[i+1] - x[i])
        y[i+1] = d*d
    end
    return y
end

function movAvg(param)
    function ewma!(x,y)
        #y = zeros(length(x))
        y[1] = x[1]
        for i in 1:length(x)-1
            y[i+1] = (1-param) * y[i] + param* x[i+1]
        end
        #@show y
        return y
    end
    return ewma!
end

#not numerically stable
function runvar(param)
    function uncertainty(x)
        y = zeros(length(x))
        S = 0
        T = 0
        for i in 1:length(x)-1
            S = S+( x[i]*x[i] )
            T = T+x[i]
            y[i+1] = (1/i) * (S-(T*T))
        end
        return y
    end
end

#taken from Art of Computer Programming
function knuthVar!(data, variance)
    #variance = zeros(length(data))
    n = 0
    mean = 0
    M2 = 0
    for x in data
        n += 1
        delta = x - mean
        mean += delta/n
        M2 = M2 + delta*(x-mean)
        variance[n] = M2/(n-1)
    end
    variance[1] = 0
    return variance
end

function rowmap!(f, M)
    n,m = size(M)
    out = zeros(n,m)
    for i=1:n
        out[i,:] = f(M[i,:], out[i,:])
    end
    #@show out
    return out
end

function id!(x, y)
    copy!(y, x)
end

#F = readdlm(STDIN, delim=" ")
extractors = [id! paddedDiffsqr! movAvg(0.5) knuthVar!]
function derivedFeatures(F, extractors)
    DV, T = size(F)
    numfeatures=length(extractors)
    M = zeros(numfeatures*DV, T)
    for i=1:numfeatures
        @show i extractors[i]
        @time M[(i-1)*DV+1:(i)*DV, :] = rowmap!(extractors[i], F)
    end
    return M
end

function test()
    F = [1 2 3 4 5; -1 -2 -3 -4 -5]
    post = derivedFeatures(F, extractors)
    @show post
end
