using Dates
using MAT
using DelimitedFiles

function get_time_in_ms()
    return convert(Dates.Millisecond, Dates.now())
end

function read_instance(m, i)
    file = matopen("instances/instance_$(m)_$(i).mat")
    A = read(file, "A")
    R = read(file, "R")
    close(file)
    n, m = size(A, 1), size(A, 2) # Gabriel uses n as the number of lines (experiments).
    s = Int(n/2)
    return A, R, m, n, s
end