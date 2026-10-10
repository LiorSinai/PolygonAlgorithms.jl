# This algorithm has been deprecated in favour of compute_graph_faces because:
# - it is more robust to improper polygons.
# - it is consistent with polygon directions.

#############################################################
##                 Segment Chaining                        ##
#############################################################

is_empty_segment(ev::SegmentEvent) = (ev.self_annotations.fill_above == false) && (ev.self_annotations.fill_below == false)
struct SegmentChainCandidate{T}
    chain_idx::Int
    match_chain_start::Bool
    segment_event::SegmentEvent{T}
end

function chain_segments(
    segments::AbstractVector{SegmentEvent{T}}
    ; atol::AbstractFloat=default_atol, check_closes::Bool=true, fuzzy_rtol::AbstractFloat=1.0
    ) where T
    # Note: if any of the regions intersect at a vertex, than this is not guaranteed to give consistent results
    # They might be joined into one region or presented as separate regions.
    # This algorithm can fail if the polygon is improper (it has lines jutting out)
    chains = Vector{SegmentEvent{T}}[]
    regions = Vector{SegmentEvent{T}}[]
    processed = Set{Segment2D{T}}()
    for event in segments
        if event.segment in processed
            continue
        end
        push!(processed, event.segment)
        candidates = SegmentChainCandidate{T}[]
        @debug("[chain_segment]: event=$event")
        @debug("[chain_segment]: candidates=$candidates")
        for (chain_idx, chain) in enumerate(chains)
            insert_matching_candidate!(candidates, chain, chain_idx, event; atol=atol)
        end
        if length(candidates) == 0 # start a new open chain
            chain = [
                SegmentEvent(event.segment, true, true, deepcopy(event.self_annotations)),
                SegmentEvent(event.segment, false, true, deepcopy(event.self_annotations))
            ]
            @debug("[chain_segment]: new chain")
            push!(chains, chain)
        elseif length(candidates) == 1 # check if it closes else append to chain
            candidate = candidates[1]
            chain = chains[candidate.chain_idx]
            if check_closes && closes_chain(chain, candidate; atol=atol)
                popat!(chains, candidate.chain_idx)
                push!(regions, chain)
                @debug("[chain_segment]: closed chain")
            else
                append_candidate!(chain, candidate)
                @debug("[chain_segment]: appended chain")
            end
        elseif length(candidates) == 2 # join two chains together
            cand1 = candidates[1]
            cand2 = candidates[2]
            chain1 = chains[cand1.chain_idx]
            chain2 = chains[cand2.chain_idx]
            append_candidate!(chain1, cand1)
            new_chain = join_chains!(chain1, chain2, cand1.match_chain_start, cand2.match_chain_start)
            chains[cand1.chain_idx] = new_chain
            @debug("[chain_segment]: combined chains")
            deleteat!(chains, cand2.chain_idx)
        else # confused
            throw("Matched segment $(candidate.segment) to more than 2 chains.")
        end
    end
    # TODO: it might be possible to close some open chains
    # - it is improper: the beginning and end is a segment(s) jutting out, so it can be closed with a segment in processing
    if check_closes
        if !isempty(chains)
            fuzzy_close!(chains, regions; atol=atol, rtol=fuzzy_rtol)
        end
        @assert isempty(chains) "There are still open chains at the end of processing all segments."
        return regions
    else
        return chains
    end
end

function insert_matching_candidate!(
    candidates::Vector{<:SegmentChainCandidate},
    chain::Vector{<:SegmentEvent},
    chain_idx::Int,
    event::SegmentEvent
    ; atol::AbstractFloat=default_atol
    )
    segment = event.segment
    is_match = false
    if is_same_point(chain[1].point, segment[1]; atol=atol)
        is_match = true
        match_chain_start = true
        match_idx = 1
    elseif is_same_point(chain[1].point, segment[2]; atol=atol)
        is_match = true
        match_chain_start = true
        match_idx = 2
    elseif is_same_point(chain[end].point, segment[1]; atol=atol)
        is_match = true
        match_chain_start = false
        match_idx = 1
    elseif is_same_point(chain[end].point, segment[2]; atol=atol)
        is_match = true
        match_chain_start = false
        match_idx = 2
    end
    if is_match
        is_start = match_idx == 1
        candidate = SegmentChainCandidate(
            chain_idx,
            match_chain_start,
            # add segment for the other point
            SegmentEvent(segment, !is_start, true, deepcopy(event.self_annotations))
        )
        push!(candidates, candidate)
    end
end

function append_candidate!(
    chain::Vector{<:SegmentEvent},
    candidate::SegmentChainCandidate
    ; atol::AbstractFloat=default_atol
    )
    if candidate.match_chain_start
        if length(chain) > 1 && 
            get_orientation(
                candidate.segment_event.point,
                chain[1].point,
                chain[2].point
                ; atol=atol
            ) == COLINEAR
            popfirst!(chain)
        end
        insert!(chain, 1, candidate.segment_event)
    else
        if length(chain) > 1 && 
            get_orientation(
                chain[end-1].point,
                chain[end].point,
                candidate.segment_event.point
                ; atol=atol
            ) == COLINEAR
            pop!(chain)
        end
        push!(chain, candidate.segment_event)
    end
end

function closes_chain(chain::Vector{<:SegmentEvent}, candidate::SegmentChainCandidate; atol::AbstractFloat=default_atol)
    candidate_point = candidate.segment_event.point
    if candidate.match_chain_start
        return is_same_point(chain[end].point, candidate_point; atol=atol)
    else
        return is_same_point(chain[1].point, candidate_point; atol=atol)
    end
end

function fuzzy_close!(
    chains::Vector{<:Vector{<:SegmentEvent}},
    regions::Vector{<:Vector{<:SegmentEvent}}
    ; atol::AbstractFloat, rtol::AbstractFloat=1.0
    )
    for idx in reverse(eachindex(chains))
        if is_fuzzy_closed(chains[idx], length(regions) + 1; atol=atol, rtol=rtol)
            push!(regions, popat!(chains, idx))
        end
    end
    regions
end

function is_fuzzy_closed(chain::Vector{<:SegmentEvent}, idx::Int; atol::AbstractFloat, rtol::AbstractFloat=1.0)
    if is_same_point(chain[1].point, chain[end].point; atol=atol)
        return true
    end
    gap = norm(chain[1].point, chain[end].point)
    path = map(event -> event.point, chain)
    gaps = norm.(path[1:(end-1)], path[2:end])
    mean_gap = sum(gaps) / length(gaps)
    if gap / mean_gap <= rtol
        @warn("Region $idx was not closed, but it has a relatively small gap and will be considered closed.
        |gap| / |mean_segment| = $gap / $mean_gap < $rtol.")
        return true
    end
    false
end

function join_chains!(chain1::Vector{<:SegmentEvent}, chain2::Vector{<:SegmentEvent}, match_chain1_start, match_chain2_start)
    # Note: with clever use of reverse! can change this to always modify chain1 in place for the same output
    if match_chain1_start && match_chain2_start
        # <--- --->
        return push!(reverse!(chain2), chain1...)
    elseif match_chain1_start && !match_chain2_start
        # <--- <----
        return push!(chain2, chain1...)
    elseif !match_chain1_start && match_chain2_start
        # ---> --->
        return push!(chain1, chain2...) 
    else # !match_chain1_start && !match_chain2_start
        # ----> <-----
        return push!(chain1, reverse!(chain2)...) 
    end
end
