function calculatedelta(model::AllometricModel, dataage, indexage::Real)
  (yname, xname, qname...) = propertynames(model.cols)
  dataindexage = deepcopy(dataage[!, [yname, xname, qname...]])
  dataindexage[!, xname] .= indexage

  mmage = modelmatrix(formula(model).rhs.terms[2:end], dataage)
  mmindexage = modelmatrix(formula(model).rhs.terms[2:end], dataindexage)

  β = coef(model)[2:end]'

  # Calculate the difference in model matrices
  Δ = mmindexage .- mmage

  # Compute the sum for each observation
  Δ = sum(β .* Δ, dims=2)[:]

  return Δ
end

function siteclassification(model::AllometricModel, dataage, indexage::Real)
  if indexage <= 0
    throw(DomainError("Index Age must be positive."))
  end
  Δ = calculatedelta(model, dataage, indexage)
  # Calculate the site classification
  ft = formula(model).lhs
  hdom = modelcols(ft, dataage)
  site = @. hdom + Δ

  if isa(ft, FunctionTerm)
    # Apply the function-specific prediction logic
    # Overwrite ŷ with the corrected values using the filtered data 'x'
    AllometricModels.predictbiascorrected!(site, [indexage], ft, model.σ²)
  end

  return round.(site, digits=1)
end

function hdomclassification(model::AllometricModel, dataage::AbstractDataFrame, indexage::Real, site::Vector{<:Real})
  if indexage <= 0
    throw(DomainError("Index Age must be positive."))
  elseif any(x -> x < 0, site)
    throw(DomainError("Site values must be positive"))
  end
  Δ = calculatedelta(model, dataage, indexage)
  (yname, xname, qname...) = propertynames(model.cols)
  sitedata = deepcopy(dataage[!, [yname, xname, qname...]])
  sitedata[!, yname] .= site
  sitedata[!, xname] .= indexage
  # Calculate the site classification
  ft = formula(model).lhs
  site = modelcols(ft, sitedata)
  hdom = @. site - Δ

  if isa(ft, FunctionTerm)
    # Apply the function-specific prediction logic
    # Overwrite ŷ with the corrected values using the filtered data 'x'
    AllometricModels.predictbiascorrected!(hdom, dataage[!, xname], ft, model.σ²)
  end

  return round.(hdom, digits=1)
end

function sitetable(model::AllometricModel, indexage::Real, hi::Real)
  if hi <= 0
    throw(DomainError("The height increment must be positive."))
  end
  # Extract property names from the fitted model's data
  (hd, age, q...) = propertynames(model.cols)
  # Calculate site classification
  site = siteclassification(model, data, indexage)
  # Get unique and sorted site classes
  sites = classcenter.(site, hi) |> unique |> sort
  # Get unique and sorted ages
  ages = model.cols[age] |> unique |> sort
  # Repeat ages for each site class
  repeatedages = repeat(ages, outer=length(sites))
  # Create repeated site classes
  repeatedsites = vcat([fill(s, length(ages)) for s in sites]...)
  # Create initial DataFrame with repeated ages and sites
  sitetable = DataFrame(age => repeatedages, hd => repeatedsites)
  # Predict the dominant heights for the site classification
  hdompredict = hdomclassification(model, sitetable, indexage, repeatedsites)
  # Insert the predicted heights into the DataFrame
  insertcols!(sitetable, :site => hdompredict, makeunique=true)
  # transform into a pivote table by site
  sitetable = unstack(sitetable, propertynames(sitetable)...) |> dropmissing!
  # Rename columns to reflect site classes
  newcolumnnames = [Symbol("S$(s)") for s in names(sitetable)[2:end]]
  rename!(sitetable, [age; newcolumnnames])
  # Create the site plot
  # uniquesites = sort(sites, rev=true)
  # sitetrace = AbstractTrace[]
  # for (i, s) in enumerate(uniquesites)
  #   idx = findall(x -> x == s, repeatedsites)
  #   push!(sitetrace, scatter(
  #     x=repeatedages[idx],
  #     y=hdompredict[idx],
  #     mode="markers+lines",
  #     name=string(s),
  #     marker=attr(
  #       opacity=0.6,
  #       linewidth=0
  #     ),
  #     line=attr(width=2)
  #   ))
  # end

  # layout = Layout(
  #   title="Site Classification Plot",
  #   xaxis=attr(title=string(age)),
  #   yaxis=attr(title=string(hd))
  # )

  # siteplot = plot(sitetrace, layout)

  # return SiteAnalysis(sitetable, siteplot)
  sitetable
end