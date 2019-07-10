a = Dict( "finame"  => "test",
            "val1"  => 1,
            "val2"  => 100)

for key in keys(a)
    println(@sprintf("%10s = %10s", key, string(a["$key"])))
end


a = [1,2,3]
b = mapreduce(identity, +, a)
b = mapreduce(x->x^2, +, a)

c = zeros(5)
fill!(c, 3)


m = 10
v = collect(1:100)
TN=10
function _f_kurt1(v, m)
    return (v - m)^2
end

z2 = map(x -> _f_kurt1(x, m), v[1:TN])

a = [1,2,3]
map(x -> x^2, a)


v = collect(1:10)
TN = 3
n = length(v)
A = Array{Float64, 2}(undef, TN, n-TN+1)
for k = TN:n
    icount = 1
    println(k)
    A[:, k-TN+1] = view(v, k-TN+1:k)
end

m = mean!(ones(n-TN+1)', A)

varm(A, m.parent', dims=1, corrected = false)

A = v

using StatsBase
c = [false, false, false, true]
a = [1,2,3,4]
mean(a, weights(c))
