function calculatedelta(model::AllometricModel, data, indexage::Real)
  ts = terms(model.formula)
  cterms = filter(t -> isa(t, ContinuousTerm), ts)

  if length(cterms) >= 3
    throw(ArgumentError("Detected multiple continuous predictor variables. Site Index requires exactly one continuous predictor."))
  end
  if !istable(data)
    throw(ArgumentError("data must be a valid table"))
  else
    cols = columntable(data)
  end

  (yname, xname, qname...) = propertynames(model.data)
  dataindexage = deepcopy(cols[[yname, xname, qname...]])
  dataindexage[xname] .= indexage
  mmage = modelmatrix(formula(model).rhs.terms[2:end], cols)
  mmindexage = modelmatrix(formula(model).rhs.terms[2:end], dataindexage)
  β = coef(model)[2:end]'
  # Calculate the difference in model matrices
  Δ = mmindexage .- mmage
  # Compute the sum for each observation
  Δ = sum(β .* Δ, dims=2)[:]

  return Δ, cols
end

function siteclassification(model::AllometricModel, data, indexage::Real)
  if indexage <= 0
    throw(DomainError("Index Age must be positive."))
  end
  Δ, cols = calculatedelta(model, data, indexage)
  # Calculate the site classification
  ft = formula(model).lhs
  hdom = modelcols(ft, cols)
  site = @. hdom + Δ

  if isa(ft, FunctionTerm)
    # Apply the function-specific prediction logic
    # Overwrite ŷ with the corrected values using the filtered data 'x'
    AllometricModels.predictbiascorrected!(site, fill(indexage, length(site)), ft, model.σ²)
  end

  return round.(site, digits=1)
end

function hdomclassification(model::AllometricModel, data, indexage::Real, site::Vector{<:Real})
  if indexage <= 0
    throw(DomainError("Index Age must be positive."))
  elseif any(x -> x < 0, site)
    throw(DomainError("Site values must be positive"))
  end
  Δ, cols = calculatedelta(model, data, indexage)
  (yname, xname, qname...) = propertynames(model.data)
  sitedata = deepcopy(cols[[yname, xname, qname...]])
  sitedata[yname] .= site
  sitedata[xname] .= indexage
  # Calculate the site classification
  ft = formula(model).lhs
  site = modelcols(ft, sitedata)
  hdom = @. site - Δ

  if isa(ft, FunctionTerm)
    # Apply the function-specific prediction logic
    # Overwrite ŷ with the corrected values using the filtered data 'x'
    AllometricModels.predictbiascorrected!(hdom, cols[xname], ft, model.σ²)
  end

  return round.(hdom, digits=1)
end

classcenter(x::Real, hi::Real) = round(x / hi) * hi + (hi / 2)

function sitetable(model::AllometricModel, indexage::Real, hi::Real)
  if hi <= 0
    throw(DomainError("The height increment must be positive."))
  end

  names_prop = propertynames(model.data)
  y_name = names_prop[1]
  age_name = names_prop[2]

  site_vals = siteclassification(model, model.data, indexage)

  sites = sort(unique(classcenter.(site_vals, hi)))
  ages = sort(unique(model.data[age_name]))
  n_ages = length(ages)

  pred_vals = []
  for name in names_prop
    if name == age_name
      push!(pred_vals, ages)
    else
      val_type = model.data[name][1]
      push!(pred_vals, fill(val_type, n_ages))
    end
  end

  pred_data = NamedTuple{names_prop}(Tuple(pred_vals))

  col_names = Symbol[age_name]
  col_data = Any[ages]

  for s in sites
    current_site_vec = fill(s, n_ages)
    hpred = hdomclassification(model, pred_data, indexage, current_site_vec)

    push!(col_names, Symbol("S$(round(s, digits=1))"))
    push!(col_data, hpred)
  end

  return NamedTuple{Tuple(col_names)}(Tuple(col_data))
end
