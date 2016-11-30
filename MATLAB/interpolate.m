function y = interpolate(table,x)
	if size(table,1) ~= 2
		error('bad table shape for interpolate')
	end
	
	points = size(table,2);
	
	if (x > table(1,points)) || (x < table(1,1))
		error('x value out of bounds for interpolation')
	end
	
	% table in order from lowest x to highest x
	% x is in first column of 2 column table
	for i = 1:points
		if table(1,i) > x
			i_hi = i;
			i_lo = i - 1;
			break
		end
	end
	
	slope = (table(2,i_hi) - table(2,i_lo)) / (table(1,i_hi) - table(1,i_lo));
	y = slope * (x - table(1,i_hi)) + table(2,i_hi);
	
end