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
