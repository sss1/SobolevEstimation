% Uses rejection sampling to generate n IID samples according to a given pdf
% supported on [0, 2*pi]

function Xs = monte_carlo_sample(pdf, upper_bound, n)

  Xs = zeros(n, 1);
  num_sampled = 0;

  while num_sampled < n

    x = rand*2*pi;
    y = rand*upper_bound;

    if y < pdf(x) % if sample should be accepted

      num_sampled = num_sampled + 1;
      Xs(num_sampled) = x;

    end

  end

end
