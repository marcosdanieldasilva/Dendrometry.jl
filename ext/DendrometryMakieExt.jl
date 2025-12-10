# module DendrometryMakieExt

using Dendrometry
using GLMakie
using Distributions
using Statistics
using LinearAlgebra
using DataFrames
using AllometricModels

Makie.@recipe(SiteCurvePlot, allometric_model, indexage, hi) do scene
  Attributes(
    colormap=:viridis,
    linewidth=2,
    alpha=0.8,
    showlabels=true,
    labelfontsize=12,
    labeloffset=5
  )
end

function Makie.plot!(plot::SiteCurvePlot)
  model = plot.allometric_model[]
  indexage = plot.indexage[]
  hi = plot.hi[]

  datant = sitetable(model, indexage, hi)
  ages = datant[1]
  ncols = length(datant)

  n_curves = ncols - 1
  cmap_sym = plot.colormap[]
  colors = cgrad(cmap_sym, n_curves, categorical=true)

  for i in 2:ncols
    yvals = datant[i]
    colsym = keys(datant)[i]
    classstr = string(colsym)[2:end]

    curve_color = colors[i-1]

    lines!(plot, ages, yvals;
      linewidth=plot.linewidth,
      alpha=plot.alpha,
      label=classstr,
      color=curve_color
    )

    if plot.showlabels[]
      text!(plot, ages[end], yvals[end];
        text=classstr,
        align=(:left, :center),
        offset=(plot.labeloffset[], 0),
        fontsize=plot.labelfontsize[],
        color=curve_color
      )
    end
  end

  ax = current_axis()
  if !isnothing(ax)
    dn = propertynames(model.data)
    ax.xlabel = string(dn[2])
    ax.ylabel = string(dn[1])
    ax.title = "Site Classification Curves (Index Age: $indexage, HI: $hi)"
  end

  return plot
end

function modelcheck(model; size=nothing, kwargs...)
  f = isnothing(size) ? Figure() : Figure(size=size)
  modelcheck!(f[1, 1], model; kwargs...)
  return f
end

function modelcheck!(pos::Union{Figure,GridPosition,GridLayout}, model;
  color=:steelblue,
  markersize=8,
  linecolor=:red,
  linewidth=2,
  alpha=0.8,
  histcolor=(:gray, 0.5),
  surfacecolor=:viridis,
  surfacealpha=0.5,
  kwargs...
)

  datanames = propertynames(model.data)
  yname = datanames[1]
  prednames = [n for n in datanames[2:end] if eltype(model.data[n]) <: Number]

  ydata = model.data[yname]
  ypred = predict(model)
  resid = residuals(model)
  ndims = length(prednames)

  d = fit(Normal, resid)
  sortedres = sort(resid)
  n = length(resid)
  p = ((1:n) .- 0.5) ./ n
  qx = quantile.(d, p)
  qy = sortedres
  q_probs = [0.25, 0.75]
  qx_q = quantile(qx, q_probs)
  qy_q = quantile(qy, q_probs)
  slope = (qy_q[2] - qy_q[1]) / (qx_q[2] - qx_q[1])
  intercept = qy_q[1] - slope * qx_q[1]
  x_line = collect(extrema(qx))
  y_line = intercept .+ slope .* x_line

  g = GridLayout(pos)

  if ndims == 1
    xname = prednames[1]
    xdata = model.data[xname]

    ax1 = Axis(g[1, 1], title="Observed vs Fitted", xlabel=string(xname), ylabel=string(yname))
    ax2 = Axis(g[1, 2], title="Residuals vs Fitted", xlabel="Fitted Values", ylabel="Residuals")
    ax3 = Axis(g[2, 1], title="Histogram of Residuals", xlabel="Residuals", ylabel="Density")
    ax4 = Axis(g[2, 2], title="Normal Q-Q Plot", xlabel="Theoretical Quantiles", ylabel="Empirical Quantiles")

    scatter!(ax1, xdata, ydata; color=color, markersize=markersize, alpha=alpha, kwargs...)
    perm = sortperm(xdata)
    lines!(ax1, xdata[perm], ypred[perm]; color=linecolor, linewidth=linewidth)

    _plotdiagnostics!(ax2, ax3, ax4, ypred, resid, qx, qy, x_line, y_line;
      color, markersize, alpha, linecolor, linewidth, histcolor, kwargs...)

  else
    stropions = string.(prednames)
    xobs = Observable(prednames[1])
    yobs = Observable(prednames[2])

    colsize!(g, 1, Relative(0.66))

    gleft = GridLayout(g[1, 1])

    Label(gleft[1, 1], "X Axis:", halign=:right)
    menux = Menu(gleft[1, 2], options=zip(stropions, stropions), default=stropions[1])
    Label(gleft[1, 3], "Y Axis:", halign=:right)
    menuy = Menu(gleft[1, 4], options=zip(stropions, stropions), default=stropions[2])

    gsliders = GridLayout(gleft[2, 1:5])

    slider_dict = Dict{Symbol,Any}()

    for (i, name) in enumerate(prednames)
      lbl = Label(gsliders[i, 1], string(name), halign=:right)

      vals = model.data[name]
      sl_range = range(extrema(vals)..., length=100)
      start_val = eltype(vals) <: Integer ? round(Int, median(vals)) : median(vals)

      sl = Slider(gsliders[i, 2], range=sl_range, startvalue=start_val)

      val_str = lift(v -> string(round(v, digits=2)), sl.value)
      val_lbl = Label(gsliders[i, 3], val_str, halign=:left)

      slider_dict[name] = sl.value

      is_active_axis = lift((x, y) -> x == name || y == name, xobs, yobs)
      is_visible = map(!, is_active_axis)

      connect!(lbl.blockscene.visible, is_visible)
      connect!(sl.blockscene.visible, is_visible)
      connect!(val_lbl.blockscene.visible, is_visible)

      on(is_active_axis; update=true) do active
        rowsize!(gsliders, i, active ? 0 : 30)
      end
    end

    ax3d = Axis3(gleft[3, 2:5],
      title="3D Fitted Model (Slice View)",
      xlabel=map(string, xobs),
      ylabel=map(string, yobs),
      zlabel=string(yname),
      azimuth=1.275 * pi
    )

    on(menux.selection) do s
      xobs[] = Symbol(s)
      reset_limits!(ax3d)
    end
    on(menuy.selection) do s
      yobs[] = Symbol(s)
      reset_limits!(ax3d)
    end

    gright = GridLayout(g[1, 2])
    ax2 = Axis(gright[1, 1], title="Residuals vs Fitted", xlabel="Fitted Values", ylabel="Residuals")
    ax3 = Axis(gright[2, 1], title="Histogram of Residuals", xlabel="Residuals", ylabel="Density")
    ax4 = Axis(gright[3, 1], title="Normal Q-Q Plot", xlabel="Theoretical Quantiles", ylabel="Empirical Quantiles")

    xdata3d = @lift(model.data[$xobs])
    ydata3d = @lift(model.data[$yobs])
    zdata3d = ydata

    scatter!(ax3d, xdata3d, ydata3d, zdata3d;
      color=color, markersize=markersize, alpha=alpha, label="Observed", kwargs...)

    all_sliders_obs = collect(values(slider_dict))

    surftuple = lift(xobs, yobs, all_sliders_obs...) do x, y, args...
      _generate_surface_grid(model, x, y, slider_dict)
    end

    surfx = @lift $surftuple[1]
    surfy = @lift $surftuple[2]
    surfz = @lift $surftuple[3]

    surfplt = surface!(ax3d, surfx, surfy, surfz;
      colormap=surfacecolor, transparency=true, alpha=surfacealpha, shading=NoShading)

    Colorbar(gleft[3, 1], surfplt, label="Prediction", width=15, height=Relative(0.7), alignmode=Outside())

    _plotdiagnostics!(ax2, ax3, ax4, ypred, resid, qx, qy, x_line, y_line;
      color, markersize, alpha, linecolor, linewidth, histcolor, kwargs...)
  end

  return g
end

function _plotdiagnostics!(ax2, ax3, ax4, ypred, resid, qx, qy, x_line, y_line;
  color, markersize, alpha, linecolor, linewidth, histcolor, kwargs...)

  scatter!(ax2, ypred, resid; color=color, markersize=markersize, alpha=alpha, kwargs...)
  hlines!(ax2, [0], color=:black, linestyle=:dash)

  hist!(ax3, resid; normalization=:pdf, color=histcolor, strokewidth=1)
  d = fit(Normal, resid)
  xrng = range(extrema(resid)..., length=100)
  lines!(ax3, xrng, pdf.(d, xrng); color=linecolor, linewidth=linewidth)

  scatter!(ax4, qx, qy; color=color, markersize=markersize, alpha=alpha, kwargs...)
  lines!(ax4, x_line, y_line; color=linecolor, linestyle=:dash, linewidth=linewidth)
end

function _generate_surface_grid(model, xvar::Symbol, yvar::Symbol, slider_values::Dict; resolution=30)
  data = DataFrame(model.data)
  xrng = range(extrema(data[!, xvar])..., length=resolution)
  yrng = range(extrema(data[!, yvar])..., length=resolution)

  baserow = Dict{Symbol,Any}()

  for n in propertynames(data)
    if n == xvar || n == yvar
      continue
    end

    if haskey(slider_values, n)
      baserow[n] = slider_values[n][]
    elseif eltype(data[!, n]) <: Number
      baserow[n] = median(data[!, n])
    else
      baserow[n] = first(data[!, n])
    end
  end

  function predict_point(x, y)
    row = copy(baserow)
    row[xvar] = x
    row[yvar] = y
    dfsingle = DataFrame([row])
    return predict(model, dfsingle)[1]
  end

  zmatrix = [predict_point(x, y) for x in xrng, y in yrng]
  return xrng, yrng, zmatrix
end

sitecurves(model, indexage, hi; kwargs...) = sitecurveplot(model, indexage, hi; kwargs...)
sitecurves!(args...; kwargs...) = sitecurveplot!(args...; kwargs...)

Makie.plottype(::AllometricModel) = SiteCurvePlot

function Makie.convert_arguments(::Type{<:SiteCurvePlot}, model::AllometricModel)
  return (model, 15.0, 2.0)
end

# Opção 1: Específica (Recomendada)
sitecurves(best, 60, 0.5)

# Opção 2: Genérica (Usa os defaults definidos no convert_arguments)
plot(best)

# Opção 3: Dashboard de Diagnóstico
modelcheck(best)

modelcheck(reg)

modelcheck(reg2)

modelcheck(reg,
  color=:cyan,          # Pontos ciano
  markersize=10,        # Pontos maiores
  linecolor=:yellow,    # Linhas de referência amarelas
  histcolor=:white      # Histograma branco
)

modelcheck(reg,
  color=:black,         # Pontos pretos
  markersize=6,         # Pontos menores e discretos
  linecolor=:black,     # Linhas pretas
  linewidth=1,          # Linhas finas
  histcolor=(:gray80, 1.0), # Cinza claro sólido
  alpha=0.6             # Leve transparência nos pontos
)

modelcheck(reg,
  surfacecolor=:magma,  # Mapa de cor 'magma' (roxo/laranja)
  surfacealpha=0.2,     # Superfície bem transparente ("fantasma")
  color=:blue           # Pontos azuis para contrastar
)

