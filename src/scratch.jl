a = Dict( "finame"  => "test",
            "val1"  => 1,
            "val2"  => 100)

for key in keys(a)
    println(@sprintf("%10s = %10s", key, string(a["$key"])))
end
