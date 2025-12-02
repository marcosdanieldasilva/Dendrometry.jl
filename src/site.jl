function _calculate_delta(model::AllometricModel, data_age, index_age::Real)
  (yname, xname, qname...) = propertynames(model.cols)
  data_index_age = deepcopy(data_age[!, [yname, xname, qname...]])
  data_index_age[!, xname] .= index_age

  mm_age = modelmatrix(formula(model).rhs.terms[2:end], data_age)
  mm_index_age = modelmatrix(formula(model).rhs.terms[2:end], data_index_age)

  β = coef(model)[2:end]'

  # Calculate the difference in model matrices
  Δ = mm_index_age .- mm_age

  # Compute the sum for each observation
  Δ = sum(β .* Δ, dims=2)[:]

  return Δ
end

function site_classification(model::AllometricModel, data_age, index_age::Real)
  if index_age <= 0
    throw(DomainError("Index Age must be positive."))
  end
  Δ = _calculate_delta(model, data_age, index_age)
  # Calculate the site classification
  ft = formula(model).lhs
  hdom = modelcols(ft, data_age)
  site = @. hdom + Δ

  if isa(ft, FunctionTerm)
    # Apply the function-specific prediction logic
    # Overwrite ŷ with the corrected values using the filtered data 'x'
    AllometricModels.predictbiascorrected!(site, [index_age], ft, model.σ²)
  end

  return round.(site, digits=1)
end

function hdom_classification(model::AllometricModel, data_age::AbstractDataFrame, index_age::Real, site::Vector{<:Real})
  if index_age <= 0
    throw(DomainError("Index Age must be positive."))
  elseif any(x -> x < 0, site)
    throw(DomainError("Site values must be positive"))
  end
  Δ = _calculate_delta(model, data_age, index_age)
  (yname, xname, qname...) = propertynames(model.cols)
  site_data = deepcopy(data_age[!, [yname, xname, qname...]])
  site_data[!, yname] .= site
  site_data[!, xname] .= index_age
  # Calculate the site classification
  ft = formula(model).lhs
  site = modelcols(ft, site_data)
  hdom = @. site - Δ

  if isa(ft, FunctionTerm)
    # Apply the function-specific prediction logic
    # Overwrite ŷ with the corrected values using the filtered data 'x'
    AllometricModels.predictbiascorrected!(hdom, data_age[!, xname], ft, model.σ²)
  end

  return round.(hdom, digits=1)
end

_class_center(x::Real, hi::Real) = round(x / hi) * hi + (hi / 2)

function site_table(model::AllometricModel, index_age::Real, hi::Real)
  if hi <= 0
    throw(DomainError("The height increment must be positive."))
  end
  # Extract property names from the fitted model's data
  (hd, age, q...) = propertynames(model.cols)
  # Calculate site classification
  site = site_classification(model, data, index_age)
  # Get unique and sorted site classes
  sites = _class_center.(site, hi) |> unique |> sort
  # Get unique and sorted ages
  ages = model.cols[age] |> unique |> sort
  # Repeat ages for each site class
  repeated_ages = repeat(ages, outer=length(sites))
  # Create repeated site classes
  repeated_sites = vcat([fill(s, length(ages)) for s in sites]...)
  # Create initial DataFrame with repeated ages and sites
  site_table = DataFrame(age => repeated_ages, hd => repeated_sites)
  # Predict the dominant heights for the site classification
  hdom_predict = hdom_classification(model, site_table, index_age, repeated_sites)
  # Insert the predicted heights into the DataFrame
  insertcols!(site_table, :site => hdom_predict, makeunique=true)
  # transform into a pivote table by site
  site_table = unstack(site_table, propertynames(site_table)...) |> dropmissing!
  # Rename columns to reflect site classes
  new_column_names = [Symbol("S_$(s)") for s in names(site_table)[2:end]]
  rename!(site_table, [age; new_column_names])
  # Create the site plot
  # unique_sites = sort(sites, rev=true)
  # site_trace = AbstractTrace[]
  # for (i, s) in enumerate(unique_sites)
  #   idx = findall(x -> x == s, repeated_sites)
  #   push!(site_trace, scatter(
  #     x=repeated_ages[idx],
  #     y=hdom_predict[idx],
  #     mode="markers+lines",
  #     name=string(s),
  #     marker=attr(
  #       opacity=0.6,
  #       line_width=0
  #     ),
  #     line=attr(width=2)
  #   ))
  # end

  # layout = Layout(
  #   title="Site Classification Plot",
  #   xaxis=attr(title=string(age)),
  #   yaxis=attr(title=string(hd))
  # )

  # site_plot = plot(site_trace, layout)

  # return SiteAnalysis(site_table, site_plot)
  site_table
end