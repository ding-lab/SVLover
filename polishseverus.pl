$input=$ARGV[0];
$output=$ARGV[1];
open(file,"$input");
open(res,">$output");
while(<file>)
{
	chomp;
	if(grep(/^chr/,$_))
	{
		@token=split(/\t/,$_);
		if(grep(/INV/,$token[4]))
		{
			print res $token[0]."\t".$token[1]."\t".$token[2]."_1"."\t".$token[3]."\t".$token[4]."\t".$token[5]."\t".$token[6]."\t".$token[7]."\t".$token[8]."\t".$token[9]."\n";
			@token1=split(/SVLEN=/,$_);
			@token2=split(/\;/,$token1[1]);
			$len=$token2[0];
			$loc=$token[1]+$len;
			# print $len."\n";
			print res $token[0]."\t".$loc."\t".$token[2]."_2"."\t".$token[3]."\t".$token[4]."\t".$token[5]."\t".$token[6]."\t".$token[7]."\t".$token[8]."\t".$token[9]."\n";
		}else{
			print res $_."\n";
		}
	}else{
		print res $_."\n";
	}
}
close(file);
close(res);
