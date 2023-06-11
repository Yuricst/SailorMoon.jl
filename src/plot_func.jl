"""
    Miscellaneous plot functions... 
"""

function plot_circle(radius, x, y, m=50)
    circle = zeros(2,m)
    thetas = LinRange(0.0, 2π, m)
    for i = 1:m
        circle[1,i] = radius*cos(thetas[i]) + x
        circle[2,i] = radius*sin(thetas[i]) + y
    end
    return circle
end