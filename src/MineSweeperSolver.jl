module MineSweeperSolver

export Board, Queue, generateBoard, cascadeOut, getEquationsUnknowns, removeBottomZeroRows,
       generateLinearSystem, updateBoard, generateAndSolveBoard, solvableBoard

typealias Board Array{Float64,2}
typealias Queue Array{(Int,Int),1}

const Empty = zeros(Int, 0, 0)

function generateBoard(N::Int, M::Int, n_mines::Int)
## generate random board of size N*M with n_mines
## add border of NaNs for simpler for-loops
## populate play board with adjacent mine counts
## use -1 for mine, non-negative integer for count of adjacent mines
    board_base = zeros(Int, N, M)
    board_play = fill(NaN, N+2, M+2)
    r = randperm(N * M)[1 : n_mines]
    board_base[r] = -1
    board_play[2:N+1, 2:M+1] = board_base
    for i in 2:N+1
        for j in 2:M+1
            if board_play[i, j] != -1
                board_play[i, j] = length(find(board_play[i-1:i+1, j-1:j+1] .== -1))
            end
        end
    end
    return board_play, board_base
end

function cascadeOut(position::(Int,Int), hidden::Board, known::Board)
## open adjacent cells and recursively cascade out to open new safe (0) cells
    n, m = position
    for i in n-1 : n+1
        for j in m-1 : m+1
            cell = hidden[i, j]
            if cell == 0 && known[i, j] != 0
                known[i, j] = cell
                cascadeOut((i,j), hidden, known)
            end
            known[i, j] = cell
        end
    end
end

function getEquationsUnknowns(hidden::Board, known::Board)
## find positions of unknown perimeter cells (unknowns)
## and positions of RHS of linear (equations)
    unknowns = (Int, Int)[]
    equations = (Int, Int)[]
    N, M = size(hidden)
    N -= 2
    M -= 2
    for i in 2:N+1
        for j in 2:M+1
            k = known[i, j]
            if k == -Inf
                n_knowns = length(find(known[i-1:i+1, j-1:j+1] .> -1))
                if n_knowns > 0
                    push!(unknowns, (i,j))
                end
            elseif k > 0
                n = length(find(known[i-1:i+1, j-1:j+1] .== -Inf))
                if n > 0
                    push!(equations, (i,j))
                end
            end
        end
    end
    return (equations, unknowns)
end

function removeBottomZeroRows(x::Matrix{Int})
    c = 0
    for i in 1:size(x,1)
        if !all( x[i,:] .== 0)
            c += 1
        else
            return x[1:c, :]
        end
    end
    return x
end

function generateLinearSystem(equations::Queue, unknowns::Queue, known::Board)
## generate matrix representation of linear system of (equations) and unknowns (unknowns)
    n_unknowns = length(unknowns)
    n_equations = length(equations)
    LS = zeros(Int, n_equations, n_unknowns + 1)
    for n in 1:n_equations
        ei, ej = equations[n]
        for m in 1:n_unknowns
            ci, cj = unknowns[m]
            if ci-1 <= ei <= ci+1 && cj-1 <= ej <= cj+1
                LS[n, m] = 1
            end
        end
        adj_mines = length(find( known[ei-1:ei+1, ej-1:ej+1] .== -1))
        LS[n, n_unknowns + 1] = known[ei, ej] - adj_mines
    end
    return removeBottomZeroRows(int(rref(LS)))    # return row reduced linear system
end

function updateBoard(rrLS::Matrix{Int}, unknowns::Queue, hidden::Board, known::Board)
## assign solutions (mine/non-mine) of linear system to unknowns
## set board cell to -1 for mine, get/set value from hidden board for non-mine, cascadeOut() for new safe cells
    for i in 1:size(rrLS, 1)
        inds_nzeros = find(rrLS[i, 1:end-1] .!= 0)
        inds_ones = find(rrLS[i, 1:end-1] .== 1)
        if inds_ones == inds_nzeros && length(inds_ones) > 0
            if rrLS[i, end] == 0                        ## non-mines
                for k in inds_ones
                    ci, cj = unknowns[k]
                    cell = hidden[ci, cj]
                    known[ci, cj] = cell
                    if cell == 0
                        cascadeOut( (ci, cj), hidden, known)
                    end
                end
            elseif rrLS[i, end] == length(inds_ones)    ## mines
                for k in inds_ones
                    ci, cj = unknowns[k]
                    known[ci, cj] = -1
                end
            end
        end
    end
end

function generateAndSolveBoard(N::Int, M::Int, n_mines::Int)
## generate hidden and known play boards, N*M with n_mines; set unknown cells to -Inf
## find safe start-cells from hidden board, if none exist, return
## randomly select a safe start-cell and cascadeOut()
## iterate: construct row-reduced linear system of equations and unknowns; update board with solutions to system
## return board if solved
    c_iter = 1
    hidden, board_base = generateBoard(N, M, n_mines)
    safe_cells = find(hidden .== 0)
    n_safe_cells = length(safe_cells)
    if n_safe_cells == 0
        println("Initialize Board: No Safe Cells")
        return Empty
    end
    cell_ind = randperm(n_safe_cells)[1]
    i, j = ind2sub(size(hidden), safe_cells[cell_ind])
    cell = hidden[i, j]
    known = fill(NaN, size(hidden))
    known[2:end-1, 2:end-1] = -Inf
    known[i, j] = cell
    println("Initialize Board:    $(N*M) unknowns")
    cascadeOut( (i,j), hidden, known)
    prev_n_cells_left = length(find( known .== -Inf ))
    print("\tCascade Out: $prev_n_cells_left unknowns")
    while true
        n_cells_left = length(find( known .== -Inf ))
        print("\n\tIteration $c_iter: ")
        equations, unknowns = getEquationsUnknowns(hidden, known)
        if length(equations) == 0 || length(unknowns) == 0
            println("$n_cells_left unknowns - Got Stuck")
            return Empty
        end
        rrLS = generateLinearSystem(equations, unknowns, known)
        updateBoard(rrLS, unknowns, hidden, known)
        n_cells_left = length(find( known .== -Inf ))
        print("$n_cells_left unknowns")
        if n_cells_left == prev_n_cells_left
            println(" - Got Stuck")
            return Empty
        elseif n_cells_left == 0
            println()
            return board_base
        end
        prev_n_cells_left = n_cells_left
        c_iter += 1
    end
end

function solvableBoard(N::Int, M::Int, n_mines::Int)
## return a solvable board of size N*M with n_mines
## values: 0 for non-mine, 1 for mine
    if N < 3 || M < 3 || n_mines < 1
        return
    end
    while true
        board = generateAndSolveBoard(N, M, n_mines)
        if board != Empty
            println("\nSolvable Board with $n_mines mines")
            return abs(board)
        end
    end
end

end # module
