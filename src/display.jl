#RSDeltaSigmaPort display facilities
#-------------------------------------------------------------------------------

const PlotOrColl = Union{EasyPlot.Plot, EasyPlot.PlotCollection}

struct ImageRepr{T}
	obj::T #Object to display as image
	width::Float64 #Width of resultant image
	AR::Float64 #Aspect ratio
end

"""`ImageRepr{T}(obj::T; width=640, AR=1)`

Wrapper object requesting an image representation of something using a specific
width and aspect ratio.
"""
ImageRepr(obj; width=640, AR=1) = ImageRepr(obj, Float64(width), Float64(AR))


#==helper functions
===============================================================================#
#Convert EasyPlot.Plot -> EasyPlot.PlotCollection
_plt2coll(plot::EasyPlot.Plot) = push!(cons(:plot_collection), plot)


#=="show" interface
===============================================================================#
function _show(io::IO, mime::MIME, r::ImageRepr{EasyPlot.PlotCollection})
	w=round(Int, r.width); h = round(Int, w/r.AR)
	set = EasyPlot.set
	opt = EasyPlot.ShowOptions(dim=set(w=w, h=h))
	b = EasyPlot.getbuilder(:image, :InspectDR)
	nativeplot = EasyPlot.build(b, r.obj)
	EasyPlot._show(io, mime, opt, nativeplot)
	return nothing
end
function _show(io::IO, mime::MIME, r::ImageRepr{EasyPlot.Plot})
	pcoll = _plt2coll(r.obj)
	return _show(io, mime, ImageRepr(pcoll, width=r.width, AR=r.AR))
end

#Only publicly provide Base.show() for MIME"image/png":
Base.show(io::IO, mime::MIME"image/png", r::ImageRepr{T}) where T<:PlotOrColl = _show(io, mime, r)

function Base.show(io::IO, mime::MIME"text/plain", r::ImageRepr{EasyPlot.PlotCollection})
	#Do not warn. Jupyter calls function even if it doesn't use the results.
	#Might it be generating text for a hover-pop-up or something?
#	@warn("Trying to display inline plot on simple text display.")
	println(io, "[INLINE PLOT]: ", r.obj.title)
	return nothing
end
function Base.show(io::IO, mime::MIME"text/plain", r::ImageRepr{EasyPlot.Plot})
	println(io, "[INLINE PLOT]: ", r.obj.title)
	return nothing
end

#=="showable" interface
===============================================================================#
#Base.showable(mime::MIME"image/png", p::ImageRepr) = true


#==User-facing interface
===============================================================================#
"""`inlinedisp(plot; AR=1, scalew=1, maxw=900)`

Display plot as inline image.

# Arguments
 - `AR`: Aspect ratio (w/h) of generated image.
 - `maxw`: Reference width of image (typ. max width of display).
 - `scalew`: Amount by which to scale image before rendering it (`imagew=maxw*scalew`)
"""
inlinedisp(plot; AR=1, scalew=1, maxw=900) =
	display(ImageRepr(plot; AR=AR, width=maxw*scalew))

_write(mime::Symbol, filepath::String, nativeplot; opt_kwargs...) =
	_write(mime, filepath, ShowOptions(; opt_kwargs...), nativeplot)


"""`saveimage(mime::Symbol, filepath::String, plot; AR=1, width=900)`

Saves plot as an image.

# Arguments
 - `AR`: Aspect ratio (w/h) of generated image.
 - `maxw`: Reference width of image (typ. max width of display).
 - `scalew`: Amount by which to scale image before rendering it (`imagew=maxw*scalew`)
"""
function saveimage(mime::Symbol, filepath::String, plot::PlotOrColl; AR=1, width=900)
	ir = ImageRepr(plot; AR=AR, width=width)
	open(filepath, "w") do io
		_show(io, EasyPlot._getmime(mime), ir)
	end
end

function displaygui(plot::PlotOrColl)
	result = EasyPlot.displaygui(:InspectDR, plot)
	sleep(0.1) #Get GUI to refresh (Ideally: `doevents()`, or something)
	#NOTE: Appears sleep needs to be > 0 to refresh.
	return result
end

#Last line
