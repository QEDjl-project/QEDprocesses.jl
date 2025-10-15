function LARGE_TESTS()
    # if LARGE_TESTS is set to 0 or 1, use that to decide
    LT = get(ENV, "LARGE_TESTS", "CI")

    if LT == "CI"
        # by default, run large tests locally but not in CI
        return !tryparse(Bool, get(ENV, "CI", "0"))
    else
        return tryparse(Bool, LT)
    end
end

"""
    _is_test_platform_active(env_vars::AbstractVector{String}, default::Bool)::Bool

# Args
- `env_vars::AbstractVector{String}`: List of the names of environment variables. The value of the
    first defined variable in the list is parsed and returned.
- `default::Bool`: If none of the variables named in `env_vars` are defined, this value is returned.

# Return

Return if platform is active or not.
"""
function _is_test_platform_active(env_vars::AbstractVector{String}, default::Bool)::Bool
    for env_var in env_vars
        if haskey(ENV, env_var)
            return tryparse(Bool, ENV[env_var])
        end
    end
    return default
end
