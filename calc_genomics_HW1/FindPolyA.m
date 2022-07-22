% function for Q6 in genomics hw 1
% finds polyA, arranged by length
function indices =FindPolyA(Genome)

    % first, need to find the indices which presumably contains a polyA
    % find out if the genome starts with A

    if strcmp(Genome(1),'a')
        start=1;
    else
        start=[];
    end

    % find starting points for potential polyA
    % if the index before is not A, and index after is A, therefore
    % index+1 is a start of a potential polyA

    ca = strfind(Genome,'ca') + 1; ca=ca';
    ga = strfind(Genome,'ga') + 1; ga=ga';
    ta = strfind(Genome,'ta') + 1; ta=ta';
    start = sort([start;ca;ga;ta]);

    % define a memory for the longest polyA checked so far
    % define the end result

    absolute_max_length = 0; indices = [];
    
    for i=1:length(start)
        counter = 1;

    % anti error: if the index is the last one in the sequence- dont go
    % over the last index of Genome and create an index error.

      if length(Genome) <= start(i) + 1
             break
      end
    
    % check how many 'a' are in the polyA
      while strcmp(Genome(start(i)+counter),'a')
             counter = counter+1;

        % anti error: if the index is the last one in the sequence- dont go
        % over the last index of Genome and create an index error.

             if length(Genome) <= start(i) + counter
                break
             end
      end

      % as we can see, 'counter' is exactly the polyA length
      % if polyA found is the new longest, take him instead of the later
      % if polyA found is same length as last one, take him as well
      % if polyA is shorter - do nothing

      if counter > absolute_max_length
              absolute_max_length = counter;
              indices = start(i);
      elseif counter == absolute_max_length
              indices = [indices,start(i)];
      end
    end
end