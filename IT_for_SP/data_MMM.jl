mutable struct DataMMM
  n::Int
  hc::Vector{Float64} # nominal cost
  d::Vector{Float64} # deviation

  C::Int
  w::Vector{Int}

  from::Vector{Int}
  to::Vector{Int}
  δ⁺::Array{Set,1}
  δ⁻::Array{Set,1}
  nV::Int
  A::Matrix{Int}
  s::Int
  t::Int

  DataMMM() = new()
end

mutable struct solution
  TTime::Int64
  Tsol::Int64
  LP::Float64
  opt::Float64
  n_sol::Int64
  n_sol_red::Int64

  solution() = new()
end

function read_dataMMM(instance::AbstractString, number::Int)
  data = DataMMM()
  f = open(string(dirname(@__FILE__),"/data/",instance))
    ln = readline(f)
    this_instance = false
    while !this_instance
      ln = readline(f)
      parse = split(ln)
      this_number = convert(Int64,float(parse[2]))
      this_instance = (this_number == number)
      if !this_instance
        for i in 1:10 ln = readline(f) end
      end
    end
    ln = readline(f)
    parse = split(ln)
    data.n = convert(Int64,float(parse[2]))
    ln = readline(f)
    #parse = split(ln)
    #data.Γ = convert(Int64,float(parse[6]))
    ln = readline(f)
    ln = readline(f)
    data.w = extract_vector(ln)
    ln = readline(f)
    parse = split(ln)
    data.C = convert(Int64,float(parse[2]))
    ln = readline(f)
    ln = readline(f)
    data.hc = extract_vector(ln)
    ln = readline(f)
    ln = readline(f)
    data.d = extract_vector(ln)
  close(f)
  return data
end

function extract_vector(line::String)
  line = replace(line, "(", "")
  line = replace(line, ")", "")
  line = replace(line, ",", " ")
  vector = Float64[]
  parse = split(line)
  for element in parse push!(vector,convert(Float64,float(element))) end
  return vector
end

function read_dataMMM_SP(instance::AbstractString)
  data = DataMMM()
  file = readdlm(string(dirname(@__FILE__),instance))
  println("********** Read instance $(file[2,1]) ****************")
  data.nV = file[4,1]
  data.n = file[6,1]
  coord = file[8:(7+data.nV),1:2]
  data.A = file[(9+data.nV):(8+data.nV+data.n),1:data.nV]
  data.s = file[10+data.nV+data.n,1]
  data.t = file[12+data.nV+data.n,1]

  println(data.n)
  data.from = Vector{Int}(undef,data.n)
  data.to = Vector{Int}(undef,data.n)

  for i in 1:data.n
    data.from[i] = findall(data.A[i,1:data.nV] .== 1)[1]
    data.to[i] = findall(data.A[i,1:data.nV] .== -1)[1]
  end

  data.δ⁺ = Array{Set,1}(undef,data.nV)
  data.δ⁻ = Array{Set,1}(undef,data.nV)
  for i in 1:data.nV
    data.δ⁺[i] = Set()
    data.δ⁻[i] = Set()
  end
  for i in 1:data.n
    push!(data.δ⁺[data.from[i]],i)
    push!(data.δ⁻[data.to[i]],i)
  end

  #println(data.δ⁺)
  #println("from $from")
  #println("to $to")
  #println(coord)

  data.d = ones(data.n)
  data.hc = ones(data.n)
  for i in 1:data.n
    data.hc[i] = norm(coord[data.from[i],1:2]-coord[data.to[i],1:2])
    data.d[i] = data.hc[i]/2
  end
  return data
end
