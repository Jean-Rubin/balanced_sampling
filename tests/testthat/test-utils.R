describe("clamp()", {
  it("leaves the value unchanged if it is between bound", {
    expect_equal(clamp(0.5, 0.1, 0.9), 0.5)
  })
  it("clamps the value to the given lower bound when it is lower than it", {
    expect_equal(clamp(0.5, 0.6, 0.9), 0.6)
  })
  it("clamps the value to the given upper bound when it is greater than it", {
    expect_equal(clamp(0.5, 0.1, 0.4), 0.4)
  })

  it("works on vector values", {
    expect_equal(clamp(c(0.1, 0.4, 0.9), 0.3, 0.8), c(0.3, 0.4, 0.8))
  })
})
