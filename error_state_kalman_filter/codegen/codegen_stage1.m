classdef codegen_stage1
	properties
		fid
	end

	methods

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
		%================================%
		% derived symbols simplification %
		%================================%

		%simplify the symbolic deriviation result
		mat = simplify(mat);

		%factor out common expressions
		[iters, optimized_mat, common_vars] = obj.optimize_deriviation(mat, '');

		%=========================%
		% optimize common factors %
		%=========================%
		disp('optimize common factors');
		further_optimize = 0;
		optimized_commons = {};
		index = 0;
		c_suffix = '_';
		complete = 0;
		while complete == 0
			%iteratively factor out common factors

			[iters, optimized_common_vars, factors_of_common] = ...
				obj.optimize_deriviation(common_vars, c_suffix);

			if iters > 1
				further_optimize = 1;
				index = index + 1;

				common_vars = factors_of_common;
				optimized_commons{index} = optimized_common_vars;

				c_suffix = append('_', c_suffix);
			end

			if iters == 1
				complete = 1;
			end

			%disp(iters);
		end

		%celldisp(optimized_commons);
		%disp(iters);
		%disp(optimized_mat);
		%disp(common_vars);

		%============================================%
		% common factors cannot be further optimized %
		%============================================%
		disp('common factors cannot be further optimized');
		if index == 0
			[row, column] = size(common_vars);
			
			for i = 1:row
				matlab_ccode = char(ccode(common_vars(i, 1)));
				my_ccode = strrep(matlab_ccode, '  t0 =', '');

				str = sprintf('float c%d =%s\n', i - 1, my_ccode);
				fprintf(obj.fid, str);
				%disp(str);
			end
			fprintf(obj.fid, "\n");
		end

		%=========================================%
		% common factors can be further optimized %
		%=========================================%
		disp('common factors can be further optimized')
		while index >= 1
			%disp(index);

			c_suffix = '';
			suffix_cnt = index;
			while suffix_cnt >= 1
				c_suffix = ['_', c_suffix];
				suffix_cnt = suffix_cnt - 1;
			end

			%save common expressions
			[row, column] = size(optimized_commons{index});

			for i = 1:row
				matlab_ccode = char(ccode(optimized_commons{index}(i, 1)));
				my_ccode = strrep(matlab_ccode, '  t0 =', '');

				str = sprintf('float c%d%s =%s\n', i - 1, c_suffix, my_ccode);
				fprintf(obj.fid, str);
				%disp(str);
			end
			fprintf(obj.fid, "\n");

			index = index - 1;
		end

		%save simplified common factor expressions
		if further_optimize == 1 
			[row, column] = size(optimized_common_vars);
			
			for i = 1:row
				matlab_ccode = char(ccode(optimized_common_vars(i, 1)));
				my_ccode = strrep(matlab_ccode, '  t0 =', '');

				str = sprintf('float c%d =%s\n', i - 1, my_ccode);
				fprintf(obj.fid, str);
				%disp(str);
			end
			fprintf(obj.fid, "\n");
		end
	
		%====================================%
		% save non-common factor expressions %
		%====================================%
		disp('save non-common factor expressions');
		[row, column] = size(optimized_mat);

		for r = 1:row
			for c = 1:column
				if isequal(optimized_mat(r, c), sym('0')) == 0
					matlab_ccode = char(ccode(optimized_mat(r, c)));
					my_ccode = strrep(matlab_ccode, '  t0 =', '');

					str = sprintf('%s(%d, %d) =%s\n', ...
						      prompt_str, r - 1, c - 1, my_ccode);

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
