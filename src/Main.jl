println("OHAYO")
println("wat")

x = 4 + 3
x

using VegaLite, VegaDatasets, Query, DataFrames
cars = dataset("cars")

cars |> @select(:Acceleration, :Name) |> collect


function foo(data, origin)
    df = data |> @filter(_.Origin == origin) |> DataFrame
    return df |> @vlplot(:point, :Acceleration, :Miles_per_Gallon)
end


p = foo(cars, "USA")
#p |> save("foo.png")

