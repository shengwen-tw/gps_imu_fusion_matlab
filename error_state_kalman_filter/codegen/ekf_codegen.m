classdef ekf_codegen
	properties
		filename
		fid
		mat_symbol_list = {}
		mat_symbol_size = 0
	end

	methods

	function ret_obj = preload_mat_symbol(obj, symbol_name_str, row, column)
		obj.mat_symbol_size = obj.mat_symbol_size + 1;
		obj.mat_symbol_list(obj.mat_symbol_size, :) = {symbol_name_str, [row, column]};
		%celldisp(obj.mat_symbol_list{obj.mat_symbol_size});
		ret_obj = obj;
	end

	function formatted_str = format_matrix_indexing(obj, raw_str)
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

					raw_str = ...
						strrep(raw_str ,orig_str, replace_str);
					%disp(raw_str);
				end
			end
		end

		formatted_str = raw_str;
	end

	function ret_obj = open_file(obj, filename)
		obj.filename = filename;
		obj.fid = fopen(filename, 'w');
		ret_obj = obj;
	end

	function close_file(obj)
		fclose(obj.fid);
		str = sprintf('codegen: %s is created', obj.filename);
		disp(str);
	end

	function add_c_comment(obj, comment)
		str = sprintf('%s\n', comment);
		fprintf(obj.fid, str);	
	end

	function [common_factor_cnt, optimized_expr, common_factors] = ...
		 optimize_deriviation(obj, expr, suffix)

		syms optimized_expr old_expr

		common_factor_cnt = 0;
		max_iter = 1000;
		old_expr = expr;
		common_factor_list = [];

		while 1
			%common variable name (string)
			common_var_name = ['c', num2str(common_factor_cnt), suffix];

			%factor out common expression
			[optimized_expr, new_common_factors] = ...
				subexpr(old_expr, common_var_name);

			common_factor_list = [common_factor_list; new_common_factors];

			%exit if no more common factor
			if isequaln(optimized_expr, old_expr) == 1
				break;
			end

			old_expr = optimized_expr;
			common_factor_cnt = common_factor_cnt + 1;

			if common_factor_cnt >= max_iter
				break;
			end
		end

		optimized_expr = optimized_expr;
		common_factors = common_factor_list;
	end

	function generate_c_code(obj, prompt_str, mat, is_symmetry)
		%=================================================%
		% factor out common factors from the given matrix %
		%=================================================%
		%disp('factor out common expressions from input matrix');
		[common_factor_cnt, optimized_mat, common_vars] = obj.optimize_deriviation(mat, '');
		%disp(optimized_mat);
		%disp(common_vars);

		no_common_factors = 0;
		if isempty(common_vars)
			%no common factors can be extracted, end of the optimization
			no_common_factors = 1;
		end

		%===================================================================%
		% further optimization: factor out common factors of common factors %
		%===================================================================%
		%disp('optimize common expressions');
		old_expr = common_vars;
		optimized_commons = {};
		optimize_depth = 0;
		c_suffix = '_';
		while no_common_factors == 0
			optimize_depth = optimize_depth + 1;

			%iteratively factor out common expression
			[common_factor_cnt, optimized_expr, new_common_factors] = ...
				obj.optimize_deriviation(old_expr, c_suffix);

			%append current optimized result to the last of the matlab cell
			optimized_commons{optimize_depth} = optimized_expr;

			%prepare for next iteration
			old_expr = new_common_factors;

			%accumulate suffix symbol '_' when new common variable is found
			if common_factor_cnt > 0
				c_suffix = append('_', c_suffix);
			end

			if common_factor_cnt == 0
				break;
			end
		end

		%disp(optimize_depth);
		%celldisp(optimized_commons);

		%====================================%
		% generate c code for common factors %
		%====================================%
		%disp('save common factors');
		while optimize_depth > 0
			%disp(optimize_depth);

			%generate suffix of the common variable name according to
			%the optimize_depth size
			c_suffix = '';
			suffix_cnt = optimize_depth - 1;
			while suffix_cnt >= 1
				c_suffix = ['_', c_suffix];
				suffix_cnt = suffix_cnt - 1;
			end

			%get list size (row majoring, column is always 1)
			[row, column] = size(optimized_commons{optimize_depth});

			%formatting and save the expression iteratively in c style
			for i = 1:row
				matlab_ccode = char(ccode(optimized_commons{optimize_depth}(i, 1)));
				my_ccode = strrep(matlab_ccode, '  t0 =', '');

				str = sprintf('float c%d%s =%s\n', i - 1, c_suffix, my_ccode);

				str = obj.format_matrix_indexing(str);

				%============================================%
				% replace divide by 2.0 with multiply by 0.5 %
				% (which may lower the accuracy)             %
				%============================================%
				str = strrep(str, '/2.0', '*0.5');

				fprintf(obj.fid, str);
				%disp(str);
			end
			fprintf(obj.fid, "\n");

			optimize_depth = optimize_depth - 1;
		end

		%============================================%
		% generate c code for non-common expressions %
		%============================================%
		%disp('save optimized expressions');
		[row, column] = size(optimized_mat);

		%formatting and save the matrix iteratively in c style
		for r = 1:row
			for c = 1:column
				%only save non-zero terms
				if isequaln(optimized_mat(r, c), sym('0')) == 0
					matlab_ccode = char(ccode(optimized_mat(r, c)));
					my_ccode = strrep(matlab_ccode, '  t0 =', '');

					%==================================================%
					% if matrix is known to be symmetry then just copy %
					% upper triangle element to lower triangle         %
					%==================================================%
					if is_symmetry == 'is_symmetry=1'
						%c > r: upper triangle
						%c = r: diagonal
						%c < r: lower triangle
						if (c - 1) >= (r - 1)
							%upper triangle & diagonal
							str = sprintf('%s(%d, %d) =%s\n', ...
								      prompt_str, r - 1, c - 1, my_ccode);

							str = obj.format_matrix_indexing(str);
						else
							%lower triangle
							str = sprintf('%s(%d, %d) = %s(%d, %d);\n', ...
								      prompt_str, r - 1, c - 1, ...
								      prompt_str, c - 1, r - 1);
						end
					else
						%format derived result
						str = sprintf('%s(%d, %d) =%s\n', ...
							      prompt_str, r - 1, c - 1, my_ccode);
						str = obj.format_matrix_indexing(str);
					end

					%============================================%
					% replace divide by 2.0 with multiply by 0.5 %
					% (which may lower the accuracy)             %
					%============================================%
					str = strrep(str, '/2.0', '*0.5');

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
