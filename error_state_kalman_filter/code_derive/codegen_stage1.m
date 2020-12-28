classdef codegen_stage1
	properties
	end

	methods
	function format_derived_result(obj, mat)
		pkg load symbolic

		mat = simplify(mat) %simplify the symbolic deriviation result

		prompt_str = inputname(1);

		[row, column] = size(mat);

		for r = 1:row
			for c = 1:column
				if mat(r, c) != 0
					str = sprintf('%s(%d, %d) = %s', prompt_str, ...
	                                               r - 1, c - 1, char(mat(r, c)));
					disp(str);
				end
			end
		end

		disp('');
	end
	end
end
