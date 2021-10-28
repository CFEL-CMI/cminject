abstract type Field
end

noise(particle, f::Field, time) = ()  # no noise term for fields required by default