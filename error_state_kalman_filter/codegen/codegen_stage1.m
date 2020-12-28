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

	function format_derived_result(obj, prompt_str, mat)
		pkg load symbolic

		%simplify the symbolic deriviation result
		mat = simplify(mat);

		[row, column] = size(mat);

		for r = 1:row
			for c = 1:column
				if mat(r, c) != 0
					str = sprintf('%s(%d, %d) = %s\n', prompt_str, ...
	                                               r - 1, c - 1, char(mat(r, c)));
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
