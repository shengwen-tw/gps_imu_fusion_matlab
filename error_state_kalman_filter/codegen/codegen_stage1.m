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
		 optimize_deriviation(obj, expr, prefix)

		syms new_sub_expr old_expr

		complete = 0;
		index = 0;
		max_iter = 100;
		complete = 0;
		old_expr = expr;
		common = [];

		while complete == 0
			%common variable name (string)
			common_var_name_str = ['c', num2str(index), prefix]

			%factor out common expression
			[new_sub_expr, new_common] = ...
				subexpr(old_expr, common_var_name_str);

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

		%factor out common expressions
		[iters, optimized_mat, common_vars] = obj.optimize_deriviation(mat, '');

		%optimize common expressions
		c_prefix = '_';
		while iters > 1
			[iters, optimized_common_vars, factors_of_common] = ...
				obj.optimize_deriviation(common_vars, c_prefix);

			common_vars = [factors_of_common; optimized_common_vars];

			if iters > 1
				c_prefix = append('_', c_prefix)
			end
			%disp(iters);
		end

		%disp(iters);
		%disp(optimized_mat);
		%disp(common_vars);

		%save common expressions
		[row, column] = size(common_vars);
		for i = 1:row
			str = sprintf('float c%d = %s;\n', ...
				      i - 1, char(common_vars(i, 1)));
			fprintf(obj.fid, str);
			%disp(str);
		end

		fprintf(obj.fid, "\n");

		%save non-common expressions
		[row, column] = size(optimized_mat);

		for r = 1:row
			for c = 1:column
				if isequal(optimized_mat(r, c), sym('0')) == 0
					matlab_ccode = char(ccode(optimized_mat(r, c)));
					my_ccode = strrep(matlab_ccode, '   t0 =', '');

					str = sprintf('%s(%d, %d) = %s\n', ...
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
