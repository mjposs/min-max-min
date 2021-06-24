# Prints all paths from 's' to 't'
function get_all_paths(data, opt, sp)

  function getAllPathsUtil(arc::Int, u::Int, t::Int, visited::BitArray{1}, path::BitArray{1}, cost::Float64)
    # Mark the current node as visited and store in path
    if cost + sp[u] < opt
      visited[u] = true
      if arc > 0
        path[arc] = true
      end
      # If current vertex is same as destination, then store current path
      if u == t
        push!(X,copy(path))
      else
        for new_arc in data.δ⁺[u]
          if !visited[data.to[new_arc]]
            cost = cost + data.hc[new_arc]
            getAllPathsUtil(new_arc, data.to[new_arc], t, visited, path, cost)
            cost = cost - data.hc[new_arc]
          end
        end
      end
      # Remove current vertex from path[] and mark it as unvisited
      if arc > 0
        path[arc] = false
      end
      visited[u] = false
      return X
    end
  end

  X = Vector{BitArray{1}}()
  path = falses(data.n)
  cost = 0.0

  # Mark all the vertices as not visited
  visited = falses(data.nV)

  # Call the recursive helper function to print all paths
  getAllPathsUtil(0, data.s, data.t, visited, path, cost)
  return X
end
