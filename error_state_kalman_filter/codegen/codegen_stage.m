classdef codegen_stage1
	properties
		fid
		mat_symbol_list = {}
		mat_symbol_size = 0
	end

	methods

	function ret_obj = preload_mat_symbol(obj, symbol_name, row, column)
		obj.mat_symbol_size = obj.mat_symbol_size + 1;
		obj.mat_symbol_list(obj.mat_symbol_size, :) = {symbol_name, [row, column]};
		%celldisp(obj.mat_symbol_list{obj.mat_symbol_size});
		ret_obj = obj;
	end

	function formatted_str = format_matrix_indexing(obj, unformatted_str)
		symbol_list_size = size(obj.mat_symbol_list);
		%celldisp(obj.mat_symbol_list)

		for i = 1:symbol_list_size(1)
			dim = obj.mat_symbol_list{i, 2};
			row = dim(1);
			column = dim(2);
			%disp(row);
			%disp(column);

			%string replacement
			mat_name = obj.mat_symbol_list{i, 1};

			%disp(mat_name);
			for r = 1:(row)
				for c = 1:(column)
					r_str = num2str(r-1);
					c_str = num2str(c-1);

					orig_str = ...
						[mat_name, r_str, c_str];
					replace_str = ...
						[mat_name, '(', r_str, ',', c_str, ')'];

					unformatted_str = ...
						strrep(unformatted_str ,orig_str, replace_str);
					%disp(unformatted_str);
				end
			end
		end

		formatted_str = unformatted_str;
	end

	function ret_obj = open_file(obj, filename)
		obj.fid = fopen(filename, 'w');
		ret_obj = obj;
	end

	function close_file(obj, filename)
		fclose(obj.fid);
	end

	function [iters, optimized_expr, common_factors] = ...
		 optimize_deriviation(obj, expr, suffix)

		syms new_sub_expr old_expr

		complete = 0;
		index = 0;
		max_iter = 100;
		complete = 0;
		old_expr = expr;
		common = [];

		while complete == 0
			%common variable name (string)
			common_var_name = ['c', num2str(index), suffix];

			%factor out common expression
			[new_sub_expr, new_common] = ...
				subexpr(old_expr, common_var_name);

			common = [common; new_common];

			%exit if no more common factor
			if isequal(new_sub_expr, old_expr) == 1
				complete = 1;
			end

			old_expr = new_sub_expr;

			index = index + 1;

			if index == (max_iter - 1)
				complete = 1;
			end
		end

		iters = index;
		optimized_expr = new_sub_expr;
		common_factors = common;
	end

	function format_derived_result(obj, prompt_str, mat)
		%simplify the symbolic deriviation result
		mat = simplify(mat);

		%=============================================================%
		% factor out common expressions and get optimized expressions %
		%=============================================================%
		%disp('factor out common expressions from input matrix');
		[iters, optimized_mat, common_vars] = obj.optimize_deriviation(mat, '');
		%disp(optimized_mat);
		%disp(common_vars);

		%====================================%
		% optimize common factor expressions %
		%====================================%
		%disp('optimize common expressions');
		old_common_vars = common_vars;
		optimized_commons = {};
		depth = 0;
		c_suffix = '_';
		complete = 0;
		while complete == 0
			depth = depth + 1;

			%iteratively factor out common expression
			[iters, optimized_common_vars, factors_of_common] = ...
				obj.optimize_deriviation(old_common_vars, c_suffix);

			%append current optimized result to the last of the matlab cell
			optimized_commons{depth} = optimized_common_vars;

			%prepare for next iteration
			old_common_vars = factors_of_common;

			%accumulate suffix symbol '_' when new common variable is found
			if iters > 1
				c_suffix = append('_', c_suffix);
			end

			%no more common factor exists, stop the loop
			if iters == 1
				complete = 1;
			end
		end

		%celldisp(optimized_commons);

		%================%
		% common factors %
		%================%
		%disp('save common factors');
		while depth >= 1
			%disp(depth);

			%generate suffix of the common variable name according to
			%the depth size
			c_suffix = '';
			suffix_cnt = depth - 1;
			while suffix_cnt >= 1
				c_suffix = ['_', c_suffix];
				suffix_cnt = suffix_cnt - 1;
			end

			%get list size (row majoring, column is always 1)
			[row, column] = size(optimized_commons{depth});

			%formatting and save the expression iteratively in c style
			for i = 1:row
				matlab_ccode = char(ccode(optimized_commons{depth}(i, 1)));
				my_ccode = strrep(matlab_ccode, '  t0 =', '');

				str = sprintf('float c%d%s =%s\n', i - 1, c_suffix, my_ccode);

				str = obj.format_matrix_indexing(str);

				fprintf(obj.fid, str);
				%disp(str);
			end
			fprintf(obj.fid, "\n");

			depth = depth - 1;
		end

		%===============================%
		% non-common factor expressions %
		%===============================%
		%disp('save optimized expressions');
		[row, column] = size(optimized_mat);

		%formatting and save the matrix iteratively in c style
		for r = 1:row
			for c = 1:column
				%only save non-zero terms
				if isequal(optimized_mat(r, c), sym('0')) == 0
					matlab_ccode = char(ccode(optimized_mat(r, c)));
					my_ccode = strrep(matlab_ccode, '  t0 =', '');

					str = sprintf('%s(%d, %d) =%s\n', ...
						      prompt_str, r - 1, c - 1, my_ccode);

					str = obj.format_matrix_indexing(str);

					fprintf(obj.fid, str);
					%disp(str);
				end
			end
		end

		fprintf(obj.fid, "\n");
		%disp('');
	end
	end
end
