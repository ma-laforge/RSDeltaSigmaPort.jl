#RSDeltaSigma: State space utilities
#-------------------------------------------------------------------------------
#SEARCH: even, odd

#==getorder
===============================================================================#
#Figure out nu, nq, and order of the system & ensure dimensions are ok
function get_nu_nq_order(u, ABCD, nlev)
	u = conv2seriesmatrix2D(u)
	nu = size(u,1)
	nq = length(nlev)
	if !(size(ABCD,2) > 2 && size(ABCD,2)==nu+size(ABCD,1))
		msg = "The ABCD argument does not have proper dimensions."
		throw(ArgumentError(msg))
	end
	order = size(ABCD,1)-nq
	return (nu, nq, order)
end


#==partitionABCD
===============================================================================#
"""`[A B C D] = partitionABCD(ABCD, m)`

Partition ABCD into A, B, C, D for an m-input state-space system.
"""
function partitionABCD(ABCD, m=nothing)
	n = 0
	if m!= nothing
		n = size(ABCD,2) - m
	else
		n = minimum(size(ABCD)) - 1
		m = size(ABCD,2) - n
	end
	r = size(ABCD,1) - n

	A = ABCD[1:n, 1:n]
	B = ABCD[1:n, n+1:n+m]
	C = ABCD[n+1:n+r, 1:n]
	D = ABCD[n+1:n+r, n+1:n+m]

	return (A, B, C, D)
end


#==stuffABCD
===============================================================================#
function _stuffABCD(::DS{:CRFB}, a,g,b,c)
	order = length(a)
	ABCD = zeros(order+1,order+2)
	(order2, isodd, Δodd) = orderinfo(order)
	Δeven = 1-Δodd

	#C=(0 0...c_n)
	#This is done as part of the construction of A, below
	#B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
	ABCD[:,order+1] = b'
	#B2 = -(a_1 a_2... a_n)
	ABCD[1:order,order+2] = -a'
	diagonal = 1:(order+2):order*(order+1)
	ABCD[diagonal] = ones(1,order)
	subdiag = diagonal[(1+Δeven):2:order] .+ 1
	ABCD[subdiag] = c[1+Δeven:2:order]
	supdiag = subdiag[(1+Δodd):length(subdiag)] .- 2
	ABCD[supdiag] = -g
	dly = (2+Δodd):2:order #row numbers of delaying integrators
	ABCD[dly,:] = ABCD[dly,:] .+ diagm(c[dly .- 1]) * ABCD[dly .- 1,:]

	return ABCD
end

#=
switch form       
    case 'CRFF'
        %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD(:,order+1) = b';
        %B2 = -(c_1 0... 0)
        ABCD(1,order+2) = -c(1);
        diagonal = 1:(order+2):order*(order+1);	% # of elements = order
        ABCD(diagonal) = ones(1,order);
        subdiag = diagonal(1:2:order-1)+1;
        ABCD(subdiag)= c(2:2:order);
        if even
            multg = 1:2:order;	% rows to have g*(following row) subtracted.
            ABCD(multg,:) = ABCD(multg,:) - diagm(g)*ABCD(multg+1,:);
        elseif order >= 3
            supdiag = diagonal(3:2:order)-1;
            ABCD(supdiag) = -g;
        end
        multc = 3:2:order;		% rows to have c*(preceding row) added.
        ABCD(multc,:) = ABCD(multc,:) + diagm(c(multc))*ABCD(multc-1,:);
        ABCD(order+1,1:2:order) = a(1:2:order);
        for i = 2:2:order
            ABCD(order+1,:) = ABCD(order+1,:) + a(i)*ABCD(i,:);
        end
        
    case 'CIFB'
        %C=(0 0...c_n)
        % This is done as part of the construction of A, below
        %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD(:,order+1) = b';
        %B2 = -(a_1 a_2... a_n)
        ABCD(1:order,order+2) = -a';
        diagonal = 1:(order+2):order*(order+1);
        ABCD(diagonal) = ones(1,order);
        subdiag = diagonal(1:order)+1;
        ABCD(subdiag)= c;
        supdiag = diagonal((2+odd):2:order)-1;
        ABCD(supdiag) = -g;
        
    case 'CIFF'
        %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD(:,order+1) = b';
        %B2 = -(c_1 0... 0)
        ABCD(1,order+2) = -c(1);
        diagonal = 1:(order+2):order*(order+1);
        ABCD(diagonal) = ones(1,order);
        subdiag = diagonal(1:order-1)+1;
        ABCD(subdiag)= c(2:end);
        %C = (a_1 a_2... a_n)
        ABCD(order+1,1:order) = a(1:order);
        supdiag = diagonal((2+odd):2:order)-1;
        ABCD(supdiag) = -g;
        
    case 'CRFBD'
        %C=(0 0...c_n)
        ABCD(order+1,order) = c(order);
        %B1 = (b_1 b_2... b_n), D=(b_n+1 0)
        ABCD(:,order+1) = b';
        %B2 = -(a_1 a_2... a_n)
        ABCD(1:order,order+2) = -a';
        diagonal = 1:(order+2):order*(order+1);
        ABCD(diagonal) = ones(1,order);
        dly = (1+odd):2:order;	% row numbers of delaying integrators
        subdiag = diagonal(dly)+1;
        ABCD(subdiag)= c(dly);
        supdiag = subdiag((1+odd):length(subdiag))-2;
        ABCD(dly,:) = ABCD(dly,:) - diagm(g)*ABCD(dly+1,:);
        if order>2
            coupl = 2+even:2:order;
            ABCD(coupl,:) = ABCD(coupl,:) + diagm(c(coupl-1))*ABCD(coupl-1,:);
        end
        
    case 'CRFFD'
        diagonal = 1:(order+2):order*(order+1);
        subdiag = diagonal(1:order-1)+1;
        supdiag = diagonal(2:order)-1;
        %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD(:,order+1) = b';
        %B2 = -(c_1 0... 0)
        ABCD(1,order+2) = -c(1);
        ABCD(diagonal) = ones(1,order);
        multc = 2:2:order;		% rows to have c*(preceding row) added.
        if order>2
            ABCD(subdiag(2:2:end)) = c(3:2:end);
        end
        if even
            ABCD(supdiag(1:2:end)) = -g;
        else
            % subtract g*(following row) from the multc rows
            ABCD(multc,:) = ABCD(multc,:) - diagm(g)*ABCD(multc+1,:);
        end
        ABCD(multc,:) = ABCD(multc,:) + diagm(c(multc))*ABCD(multc-1,:);
        % C
        ABCD(order+1,2:2:order) = a(2:2:order);
        for i = 1:2:order
            ABCD(order+1,:) = ABCD(order+1,:) + a(i)*ABCD(i,:);
        end
        % The above gives y(n+1); need to add a delay to get y(n).
        % Do this by augmenting the states. Note: this means that
        % the apparent order of the NTF is one higher than it acually is.
        [A B C D] = partitionABCD(ABCD,2);
        A = [A zeros(order,1); C 0];
        B = [B; D];
        C = [zeros(1,order) 1];
        D = [0 0];
        ABCD = [A B; C D];
        
    case 'PFF'
        %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD(:,order+1) = b';
        odd_1 = odd;		% !! Bold assumption !!
        odd_2 = 0;			% !! Bold assumption !!
        gc = g.* c(1+odd_1:2:end);
        theta = acos(1-gc/2);
        if odd_1
            theta0 = 0;
        else
            theta0 = theta(1);
        end
        order_2 = 2*length(find(abs(theta-theta0)>0.5));
        order_1 = order - order_2;
        %B2 = -(c_1 0...0 c_n 0...0)
        ABCD(1,order+2) = -c(1);
        ABCD(order_1+1,order+2) = -c(order_1+1);
        diagonal = 1:(order+2):order*(order+1);	% # of elements = order
        ABCD(diagonal) = ones(1,order);
        i = [1:2:order_1-1 order-order_2+1:2:order]
        subdiag = diagonal(i)+1;
        ABCD(subdiag)= c(i+1);
        if odd_1
            if order_1 >= 3
                supdiag = diagonal(3:2:order_1)-1;
                ABCD(supdiag) = -g(1:(order_1-1)/2);
            end
        else
            multg = 1:2:order_1;	% rows to have g*(following row) subtracted.
            ABCD(multg,:) = ABCD(multg,:) - diagm(g(1:order_1/2))*ABCD(multg+1,:);
        end
        if odd_2
            if order_2 >= 3
                supdiag = diagonal(order_1+2:2:order)-1;
                ABCD(supdiag) = -g(1:(order_1-1)/2);
            end
        else
            multg = order_1+1:2:order; % rows to have g*(following row) subtracted.
            gg = g((order_1-odd_1)/2+1:end);
            ABCD(multg,:) = ABCD(multg,:) - diagm(gg)*ABCD(multg+1,:);
        end
        % Rows to have c*(preceding row) added.
        multc = [3:2:order_1 order_1+3:2:order];
        ABCD(multc,:) = ABCD(multc,:) + diagm(c(multc))*ABCD(multc-1,:);
        % C portion of ABCD
        i = [1:2:order_1 order_1+1:2:order];
        ABCD(order+1,i) = a(i);
        for i = [2:2:order_1 order_1+2:2:order]
            ABCD(order+1,:) = ABCD(order+1,:) + a(i)*ABCD(i,:);
        end
        
    case 'Stratos'
        % code copied from case 'CIFF':
        %B1 = (b_1 b_2... b_n), D=(b_(n+1) 0)
        ABCD(:,order+1) = b';
        %B2 = -(c_1 0... 0)
        ABCD(1,order+2) = -c(1);
        diagonal = 1:(order+2):order*(order+1);
        ABCD(diagonal) = ones(1,order);
        subdiag = diagonal(1:order-1)+1;
        ABCD(subdiag)= c(2:end);
        % code based on case 'CRFF':
        multg = 1+odd:2:order-1;	% rows to have g*(following row) subtracted.
        ABCD(multg,:) = ABCD(multg,:) - diagm(g)*ABCD(multg+1,:);
        % code copied from case 'CIFF':
        %C = (a_1 a_2... a_n)
        ABCD(order+1,1:order) = a(1:order);
        
    case 'DSFB'
        I = eye(order);
        A0 = I;
        A1 = zeros(order,order);
        indices = 2:order+1:(order-1)*order;    % subdiagonal
        A1(indices) = c(1:order-1);
        indices = order+1 : 2*(order+1) : order*order-1; %supdiagonal
        A1(indices) = -g;
        C0 = I;
        C1 = I;
        D0 = zeros(order,2);
        if odd
            C1(order,order) = 0;
        else
            A1(order-1:order,order-1:order) = [ 0 -g(end)/2; c(end-1) 0];
            C0(order,order) = 0;
            D0(order-1,2) = g(end)/2;
        end
        M = inv(I - A1*C1);
        A = M * (A0 + A1*C0 );
        B = M * ( [b(1:order)' -a'] + A1*D0 );
        C = [zeros(1,order-1) c(order)];
        D = [b(order+1) 0];
        ABCD = [A B; C D];

    otherwise
        fprintf(1,'%s error: Form %s is not yet supported.\n',mfilename,form)
end
=#

"""`ABCD = stuffABCD(a,g,b,c,form='CRFB')`

Compute the ABCD matrix for the specified structure.
See realizeNTF.m for a list of supported structures.
`mapABCD` is the inverse function.
"""
function stuffABCD(a,g,b,c, form::Symbol=:CRFB)
	order = length(a)
	if length(b)==1
		b = [b zeros(1,order)]
	end

	return _stuffABCD(DS(form), a, g, b, c)
end


#==mapABCD
===============================================================================#
#Common implementation for CRFB, CIFB, & CRFBD:
function _mapABCD_CRFB_CIFB_CRFBD(ABCD, subdiag, supdiag, order::Int, form::Symbol)
	order2, isodd, Δodd = orderinfo(order)
	Δeven = 1 - Δodd
	c = ABCD[subdiag]
	g = -ABCD[supdiag]
	if :CRFB == form
		dly = (2+Δodd):2:order #row numbers of delaying integrators.
		ABCD[dly,:] = ABCD[dly,:] .- diagm(c[dly .- 1])*ABCD[dly .- 1,:]
	elseif :CRFBD == form
		dly = (1+Δodd):2:order #row numbers of delaying integrators.
		ABCD[dly,:] = ABCD[dly,:] .+ diagm(g)*ABCD[dly .+ 1,:]
		if order>2
			coupl = 2+Δeven:2:order
			ABCD[coupl,:] = ABCD[coupl,:] .- diagm(c[coupl .- 1])*ABCD[coupl .- 1,:]
		end
	end
	a = -ABCD[1:order,order+2]'
	b = ABCD[:,order+1]'
	return (a,g,b,c)
end

_mapABCD(::DS{:CRFB}, ABCD, diagonal, subdiag, supdiag, order::Int) =
	_mapABCD_CRFB_CIFB_CRFBD(ABCD, subdiag, supdiag, order, :CRFB)
_mapABCD(::DS{:CIFB}, ABCD, diagonal, subdiag, supdiag, order::Int) =
	_mapABCD_CRFB_CIFB_CRFBD(ABCD, subdiag, supdiag, order, :CIFB)
_mapABCD(::DS{:CRFBD}, ABCD, diagonal, subdiag, supdiag, order::Int) =
	_mapABCD_CRFB_CIFB_CRFBD(ABCD, subdiag, supdiag, order, :CRFBD)

#=
case 'CRFF'
    c = [-ABCD(1,order+2) ABCD(subdiag(1:end-1))];
    g = -ABCD(supdiag);
    if even
        multg = 1:2:order;	%Rows to have g*(following row) added.
        ABCD(multg,:) = ABCD(multg,:)+diag(g)*ABCD(multg+1,:);
    end
    multc=3:2:order;	%Rows to have c*(preceding row) subtracted.
    ABCD(multc,:)=ABCD(multc,:)-diag(c(multc))*ABCD(multc-1,:);
    a(2:2:order)=ABCD(order+1,2:2:order);
    for i=2:2:order                   %Recover a coeff.
        ABCD(order+1,:)=ABCD(order+1,:)-a(i)*ABCD(i,:);
    end
    a(1:2:order)=ABCD(order+1,1:2:order);
    b = ABCD(:,order+1)';

case 'CRFFD'
    % CRFFD has an extra order. Correct the common variables
    order = order-1;
	order2, isodd, Δodd = orderinfo(order)
    odd = rem(order,2); even = ~odd;
    diagonal = diagonal(1:order);
    subdiag = diagonal(1:order-1)+1;
    supdiag = diagonal(2+Δodd:2:order)-1;
    g = -ABCD(supdiag);
    c = [-ABCD(1,order+3) ABCD(subdiag)];
    a = zeros(1,order);
    for i=1:2:order                   %Recover the a coefficients
		a(i) = ABCD(order+1,i);
		ABCD(order+1,:) = ABCD(order+1,:) - a(i)*ABCD(i,:);
    end
    a(2:2:order) = ABCD(order+1,2:2:order);
    b = ABCD(1:order+1,order+2)';
    for i=2:2:order                   %Recover the b coefficients
		b(i) = b(i) - c(i) *b(i-1);
		if isodd
		    b(i) = b(i) + g(i/2) *b(i+1);
		end
    end
    yscale =  ABCD(order+2,order+1);
    a = a*yscale;
    b(end) = b(end)*yscale;

case {'CIFF','Stratos'}
    a = ABCD(order+1,1:order);
    c = [-ABCD(1,order+2) ABCD(subdiag(1:end-1))];
    g = -ABCD(supdiag);
    b = ABCD(:,order+1)';

otherwise
    fprintf(1,'%s error: Form %s is not yet supported.\n',mfilename,form)
end
=#

"""`(a,g,b,c) = mapABCD(ABCD, form=:CRFB)`

Compute the coefficients for the specified structure.
See `realizeNTF` for a list of supported structures.
`stuffABCD` is the inverse function.
"""
function mapABCD(ABCD, form::Symbol=:CRFB)
	order=size(ABCD,1)-1
	order2, isodd, Δodd = orderinfo(order)
	diagonal = 1:(order+2):order*(order+1)
	subdiag = diagonal[1:order] .+ 1
	supdiag = diagonal[2+Δodd:2:order] .- 1

	#Algorithm below modifies ABCDs copy so as not to modify original:
	ABCD = deepcopy(ABCD)

	#Do the mapping.
	return _mapABCD(DS(form), ABCD, diagonal, subdiag, supdiag, order)
end


#==scaleABCD
===============================================================================#
"""`(ABCDs,umax,S)=scaleABCD(ABCD; nlev=2,f=0,xlim=1,ymax=nlev+5,umax,N=10^5,N0=10)`

Scale the loop filter of a general delta-sigma modulator for dynamic range.

 - ABCD: The state-space description of the loop filter.
 - nlev: The number of levels in the quantizer.
 - xlim: A vector or scalar specifying the limit for each state variable.
 - ymax: The stability threshold. Inputs that yield quantizer inputs above ymax
   are considered to be beyond the stable range of the modulator.
 - umax: The maximum allowable input amplitude. umax is calculated if it
   is not supplied.
 - ABCDs: The state-space description of the scaled loop filter.
 - S: The diagonal scaling matrix, ABCDs = [S*A*Sinv S*B; C*Sinv D];
	xs = S*x;
"""
function scaleABCD(ABCD; nlev::Int=2, f=0, xlim=1, ymax=NaN, umax=NaN, N::Int=10^5, N0=10)
	N_sim=N
	order = size(ABCD,1)-1
	if isnan(ymax)
		ymax=nlev+5
	end
	if length(xlim)==1
		xlim = xlim * ones(1,order) #Convert scalar xlim to a vector
	end
	quadrature = !isreal(ABCD)

	#Make this function repeatable
	#rng_state = rng; rng('default') #Original code
	Random.seed!(0)

	#Envelope for smooth start-up 
	raised_cosine = 0.5*(1 .- cos.(pi/N0 * (0:N0-1)))'
	if isnan(umax)
		#Simulate the modulator with DC or sine wave inputs to detect its stable
		#input range. 
		#First get a rough estimate of umax.
		ulist =(0.1:0.1:1.0) * (nlev-1)
		umax = nlev-1
		N = 10^3
		u0 = [ exp.(2π*j*f*(-N0:-1))' .* raised_cosine exp.(2π*j*f*(0:N-1))' ] .+ 0.01*[1 j]*randn(2, N+N0)
		if !quadrature; u0 = real.(u0); end
		for u in ulist
			if !quadrature
				(v, x, xmax, y) = simulateDSM(u*u0,ABCD, nlev=nlev, trackmax=true)
			else
				(v, x, xmax, y) = simulateQDSM(u*u0,ABCD, nlev=nlev, trackmax=true)
			end
			if maximum(abs.(y)) > ymax
				umax = u #umax is the smallest input found which causes 'instability'
				break
			end
		end
		if umax == ulist[1]
			umaxstr = @sprintf("%.1f", umax)
			msg = "Modulator is unstable even with an input amplitude of $umaxstr."
			throw(msg)
		end
	end
	#More detailed simulation
	N = N_sim
	u0 = [ exp.(2π*j*f*(-N0:-1))' .* raised_cosine exp.(2π*j*f*(0:N-1))' ] + 0.01*[1 j]*randn(2, N+N0)
	if !quadrature; u0 = real.(u0); end
	maxima = zeros(1,order) .- 1
	ulist = range(0.7*umax, stop=umax, length=10)

	for u in ulist
		if !quadrature
			(v, x, xmax, y) = simulateDSM(u*u0,ABCD, nlev=nlev, trackmax=true)
		else
			(v, x, xmax, y) = simulateQDSM(u*u0,ABCD, nlev=nlev, trackmax=true)
		end
		if maximum(abs.(y)) > ymax
			break
		end
		#We need to do this at least once.
		umax = u #umax is the largest input which keeps |y| < ymax
		#maxima = maximum([maxima; xmax'], dims=1)
		maxima = max.(maxima, xmax')
	end

	#Scale the modulator so that all states are at most xlim.
	scale = maxima ./ xlim
	scale = scale[:] #Convert to simple vector
	S = diagm(1 ./ scale); Sinv = diagm(scale) #xs = S * x
	(A, B, C, D) = partitionABCD(ABCD)
	ABCDs = [S*A*Sinv S*B; C*Sinv D]

	Random.seed!() #original code restores state: rng(rng_state)
	return (ABCDs, umax, S)
end


#Last line
