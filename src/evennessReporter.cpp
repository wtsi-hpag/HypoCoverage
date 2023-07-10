#include "JSL/JSL.h"

double TailTruncateFactor = 3;

struct DataFile
{
	bool GoodFile;
	std::string FileName;
	std::string Name;
	std::vector<int> Coverage;
	std::vector<double> Vals;
	double Score;
	double Mean;
	double Peak;
	double Median;
	double Variance;
	double ScoreDelta;
	std::vector<double> TruncationPoints;
	std::vector<double> TruncationValues;
	void Plot(JSL::gnuplot & gp,int i)
	{
		if (i < 12)
		{
			gp.Plot(Coverage,Vals,JSL::LineProperties::Legend(Name),JSL::LineProperties::PenSize(3));
		}
		else
		{
			gp.Plot(Coverage,Vals,JSL::LineProperties::Legend(Name),JSL::LineProperties::PenSize(3),JSL::LineProperties::PenType(JSL::DashDot));
		}
		// gp.Plot(TruncationPoints, TruncationValues,JSL::LineProperties::Legend(Name),JSL::LineProperties::PenSize(3));
	}

	void TruncationPlot(JSL::gnuplot & gp, int i)
	{
		if (i < 12)
		{
			gp.Plot(TruncationPoints, TruncationValues,JSL::LineProperties::Legend(Name),JSL::LineProperties::PenSize(3));
		}
		else
		{
			gp.Plot(TruncationPoints, TruncationValues,JSL::LineProperties::Legend(Name),JSL::LineProperties::PenSize(3),JSL::LineProperties::PenType(JSL::DashDot));
		}
		
	}

	void ComputeScore()
	{
		Mean = 0;
		double mV = 0;
		Peak = 0;
		bool medianFound = false;
		double run = 0;
		for (int i = 0; i < Coverage.size(); ++i)
		{
			Mean += Vals[i] * Coverage[i];
			if (Vals[i] > mV)
			{
				mV = Vals[i];
				Peak = Coverage[i];
			}
			run += Vals[i];
			if (run > 0.5 && !medianFound)
			{
				medianFound = true;
				if (i == 0)
				{
					Median = Coverage[i];
				}
				else
				{
					Median = Coverage[i] - (run-0.5)/(Vals[i]) * (Coverage[i] - Coverage[i-1]);
				}
			}
		}

		Variance = 0;
		for (int i =0; i < Coverage.size(); ++i)
		{
			Variance += Vals[i] * pow(Coverage[i] - Mean,2);
		}
		Score = Variance/Mean - 1.0;

		int N = 20;
		TruncationPoints = JSL::Vector::linspace(2,10,N);
		TruncationValues.resize(N);
		ScoreDelta = (Score-TruncatedScore(3))/Score;
		for (int i = 0; i < N; ++i)
		{
			TruncationValues[i]= TruncatedScore(TruncationPoints[i]);
			
		}
	}

	double TruncatedScore(double truncFactor)
	{
		double kMax = Median*truncFactor;
		double renorm = 0;
		for (int i =0; i < Coverage.size(); ++i)
		{
			if (Coverage[i] < kMax)
			{
				renorm += Vals[i];
			}
		}
		
		double truncMean = 0;
		double truncVar = 0;
		for (int i = 0; i < Coverage.size(); ++i)
		{
			Vals[i] /= renorm;
			if (Coverage[i] < kMax){
				double p = Vals[i];
				truncMean += p * Coverage[i];
				truncVar += p * Coverage[i] * Coverage[i];
			}
			
		}
		truncVar -= truncMean * truncMean;

		return (truncVar - truncMean)/truncMean;

	}

	DataFile(std::string fileLoc, const std::vector<std::vector<std::string>> & nameRefs)
	{
		// DataFile out;
		auto q = JSL::split(fileLoc,'/');

		// fileLoc = 
		FileName = q[q.size()-1];
		Name = FileName;
		if (nameRefs.size() > 0)
		{
			int q = JSL::FindXInY(fileLoc,nameRefs[0]);
			if (q > -1)
			{
				Name = nameRefs[1][q];
				std::cout << "Assigning ref " << Name << " to " << fileLoc << std::endl;
			}
			else
			{
				Name = FileName;
				std::cout << "Could not find " << fileLoc << " in name reg" << std::endl;
			}
		}
		
		GoodFile = false;
		double runSum = 0;
		forLineVectorIn(fileLoc,' ',
			if (FILE_LINE_VECTOR.size() != 2)
			{
				return;
			}
			int k = std::stoi(FILE_LINE_VECTOR[0]);
			double v = std::stod(FILE_LINE_VECTOR[1]);
			Coverage.push_back(k);
			Vals.push_back(v);
			runSum += v;
		);
		for (int i = 0; i < Vals.size(); ++i)
		{
			Vals[i] /= runSum;
		}
		GoodFile = true;
		

		ComputeScore();
	}
};
std::string getReportString(std::vector<DataFile> & files)
{
	// std::string s = "";
	std::ostringstream os;
	for (int i = 0; i < files.size(); ++i)
	{
		std::string alias = "-";
		if (files[i].FileName != files[i].Name)
		{
			alias = files[i].Name;
		}
		std::string gap = " & \\small ";
		os << "\\footnotesize " << files[i].FileName <<  gap << alias << gap;
		os << files[i].Mean << gap << files[i].Median << gap<< files[i].Variance << gap << files[i].Score << gap << files[i].ScoreDelta<<"\n\n\\tabularnewline\\hline\n";
	}
	return os.str();
}

std::string getTableString(std::vector<DataFile> & files)
{
	// std::string s = "";
	std::ostringstream os;
	os << "File,Alias,Mean Coverage, Median Coverage, Coverage Variance, Unevenness, Tail Bias\n";
	for (int i = 0; i < files.size(); ++i)
	{
		std::string alias = "-";
		if (files[i].FileName != files[i].Name)
		{
			alias = files[i].Name;
		}
		std::string gap = ",";
		os << files[i].FileName <<  gap << alias << gap;
		os << files[i].Mean << gap << files[i].Median << gap<< files[i].Variance << gap << files[i].Score << gap << files[i].ScoreDelta<<"\n";
	}
	return os.str();
}

void JFG_Latex(std::string filename, std::string compilename, bool openfile, bool quietmode)//a C++ implementation of my jfglatex bash script, to prevent carting it around
{
	std::string texTail = ".tex";
	std::string::size_type i = filename.find(texTail);
	if (i != std::string::npos)
	{
		filename.erase(i,texTail.length());
	}
	std::string cmd = "pdflatex " + filename;
	if (quietmode)
	{
		cmd += " > /dev/null"; //pipes to null to prevent printing
	}
	cmd += "; ";
	// cmd += "biblatex" //do not need biblatex for this implementation
	cmd += cmd; //double up compilation to get references working
	if (compilename.find(".pdf") == std::string::npos)
	{
		compilename += ".pdf";
	}
	
	cmd += "mv " + filename + ".pdf " + compilename + ";";
	if (JSL::locationExists("_build"))
	{
		cmd = "mv _build/* ./; " + cmd;
		// system(mvcmd.c_str());
	}
	else
	{
		JSL::mkdir("_build");
	}

	if (openfile)
	{
		cmd += " open " + compilename + "; ";
	}

	std::vector<std::string> buildArray={"*.aux", "*.log", "*.toc", "*.mtc*", "*.blg", "*.bbl", "*.out", "*.maf", "*.fdb_latexmk", "*.fls", "*-blx.bib", "*xml", ".snm", "*.nav"};

	for (int i = 0; i < buildArray.size(); ++i)
	{
		cmd += "[ -f " + buildArray[i] + " ] && mv " + buildArray[i] + " _build/;";
		// system(cmd.c_str());
	}
	system(cmd.c_str());
}

void PrintReport(std::string outputFile,std::string title,bool sortActive,std::vector<DataFile> & files, bool verbose, bool openreport)
{
	if (sortActive)
	{
		std::sort(files.begin(),files.end(),[](const DataFile & lhs, const DataFile & rhs){return lhs.Score < rhs.Score;});
	}
	if (title == "__none__")
	{
		title = "Evenness Report";
	}
	std::vector<std::string> lines = {
		"\\documentclass[10pt]{article}",
		"\\title{" + title + "}",
		"\\usepackage{amsmath}",
		"\\usepackage{graphicx}"
		"\\usepackage[left=0.8in,right=0.8in,top=0.9in,bottom=0.9in]{geometry}",
		"\\begin{document}",
			"\\maketitle",
			"\\section{Evenness Report}",
			"\\def\\s{1.8}",
			"\\def\\header{\\small \\bf \\centering}",
			"\\newcommand\\hedbox[1]{\\parbox[t]{\\s cm}{\\header #1}}"
			"\\begin{center}\\begin{tabular}{|p{2cm}|c|c|c|c|c|c|}\\hline",
				"\\header File	& \\header Alias  & \\hedbox{Mean Coverage} &\\hedbox{Median Coverage} & \\hedbox{Coverage Variance} & \\hedbox{Unevenness $\\mathcal{U}$} & \\hedbox{Tail Bias}, $ \\mathcal{T}$",
				"\\tabularnewline\\hline\\hline",
				getReportString(files),
			" \\end{tabular}\\end{center}",
			"\\begin{figure}[t] \\includegraphics[width=\\linewidth,keepaspectratio=true]{" + outputFile + "_plot.eps}\\caption{A plot of the base coverage profiles of the provided datasets}\\end{figure}",
			"\\section{Report Explained}",
			"The theoretical distribution of a coverage profile is a Poisson distribution,",
			"\\begin{equation} p(k|\\lambda) = \\frac{\\lambda^k \\exp(-\\lambda)}{k!}\\end{equation}",
			"However, it is well known that there many various processes which alter this theoretical distribution -- the most obvious being a bias towards the high-coverage tail due to",
			"over-alignment to repetitive regions. Of more concern, however, is that even near the peak of the distribution, it is frequently observed that the distribution is `fatter' than a Poisson distribution allows. This is highly indicative of multiple Poisson processes, each contributing with a different value of $\\lambda$.",
			"\\par",
			"We suggest that the observed distribution is a marginalisation of a Poisson distribution over an underlying `sampling function', $f(\\lambda)$, such that:",
			"\\begin{equation} p(k|f) = \\int_0^\\infty f(\\lambda) p(k|\\lambda) \\mathrm d \\lambda \\end{equation}",
			"If $f(\\lambda)$ is a `fat' function, or has multiple peaks, the result is a highly biased dataset, with the coverage varying in statisically unsound ways, which might lead to statistical biases further down the pipeline. We should therefore prefer $f$ to be a tightly peaked function."
			"\\par",
			"A full analysis of the behaviour of $f$ is difficult (and a work in progress at the time of writing), but it suffices to note that several properties of $f$ can be inferred simply from the observed values of $p(k|f)$. The following relationships follow automaticall from our definitions of $p(k|f)$, which we set equal to the observed distribution of coverage profiles:",
			"\\begin{align}",
				"\\langle k \\rangle & = \\langle f \\rangle \\\\",
				"\\text{Var}(f) & = \\text{Var}(k) - \\langle k \\rangle",
			"\\end{align}",
			"That is, simply by studying the mean and variance of the observed coverage profile, we can infer some statistical properties of $f$. ",
			"We define an \\textit{unevenness metric}, $\\mathcal{U}$, such that:",
			"\\begin{equation} \\mathcal{U} = \\frac{\\text{Var}(f)}{\\langle f \\rangle} = \\frac{\\text{Var}(k)}{\\langle k \\rangle} - 1 \\end{equation}",
			"This is a standard statistical dispersion metric; the coefficient of variation.,"
			"\\par",
			"In short:",
			"\\begin{itemize}",
			"\\item When $\\mathcal{U}$ is \\textbf{small}, the coverage is more Poisson-like",
			"\\item When $\\mathcal{U}$ is \\textbf{large}, the coverge is biased in some way",
			"\\end{itemize}",
			"Libraries with smaller values of $\\mathcal{U}$ should therefore be understood to have a smaller sampling bais."
			"\\subsection{Tail Weighting}"
				"One potential cause of a high unevenness score is that the data possesses a long high-coverage tail due to, i.e., over-alignment to repetitive regions. Whilst this is a statistical bias, it is an expected one, and might not be overly concerning. We therefore also introduce a metric by which can judge how large the contribution from the tail might be.",
				"\\par",
				"We compute the median coverage of the dataset, $\\bar{k}$, and define an upper limit on the dataset as $k_\\text{max} = F \\bar{k}$. For each value of $F$, we remove all values of $k$ larger than $k_\\text{max}$ (in essence, cutting off the tail at different points) and recompute $\\mathcal{U}$(F).",
				"\\par The `tail weighting' is then a measure of how much $\\mathcal{U}$ changes as the data is truncated:"
				"\\begin{equation} \\mathcal{T} = \\frac{\\mathcal{U}(\\infty) - \\mathcal{U}(3)}{\\mathcal{U}(\\infty)} \\end{equation}",
				"This value can then be interpreted as the fraction of the uneveness which arises due to the high-coverage tail. The value of 3 was chosen on the basis that a Poisson distribution with mean $>10$ should have less than 3.3\\% of its cumulative density above 3 times the median, a sufficiently small number. Figure \\ref{F:TailWeighting} shows how $\\mathcal{U}$ varies with different values of $F$.",
				"\\begin{figure}[t] \\includegraphics[width=\\linewidth,keepaspectratio=true]{" + outputFile + "_plot_variance.eps}\\caption{A plot of the value of $\\mathcal{U}$ as a function of the truncation - removing data above $F \\bar{k}$, demonstrating the effect of the high-coverage tail.}\\label{F:TailWeighting}\\end{figure}",
			"\\subsection{Data Collation}"
				"The mechanism of this metric is to detect `composite coverage', when there are multiple sources of peaks in the coverage, which are usually indicative of some hidden process in the sequencing pipeline.",
				"\\par However, in some cases this composite nature might be due to a deliberate concatenation of multiple datasets. If a 30x and a 60x sequencing effort are collated, the result is a bimodal coverage profile with a mean of 45x coverage. Even if both of these sequencing efforts were themselves perfectly even and statistically valid, the uneveness metric would flag this as a highly uneven profile (with a score of 10.0 in this case).",
				"\\par This mechanism is therefore unsuited for studying composite datasets -- it can, however, be applied to the individual datasets themselves, and if these are even, then the dataset is sound. Future work in this field will focus on deconvolving the full functional form of the sampling function $f$, and hence be suitable for composite datasets."
		"\\end{document}"
	};

	JSL::initialiseFile("tmp.tex");
	JSL::writeVectorToFile("tmp.tex",lines,"\n",false);
	JFG_Latex("tmp.tex",outputFile,openreport,!verbose);
	std::string cmd = "rm tmp.tex";
	system(cmd.c_str());
}




std::vector<std::vector<std::string>> getCrossRefs(std::string file)
{
	std::vector<std::vector<std::string>> out;
	if (file != "__none__")
	{
		out.resize(3);
		forLineVectorIn(file,' ',
			
			out[0].push_back(FILE_LINE_VECTOR[0]);
			if (FILE_LINE_VECTOR.size() > 1)
			{
				std::string s = FILE_LINE_VECTOR[1];
				for (int i = 2; i < FILE_LINE_VECTOR.size(); ++i)
				{
					s += FILE_LINE_VECTOR[i] + " ";
				}
				out[1].push_back(s);
			}
			else
			{
				auto q = JSL::split(FILE_LINE_VECTOR[0],'/');

				// fileLoc = 
				out[1].push_back(q[q.size()-1]);
			}
		);
	}
	return out;
}
void optPrint(std::ostringstream & os, std::string opt, bool keyval, std::string desc)
{
	std::string kStr = "toggle";
	if (keyval)
		kStr = "keyval";
	os << "-" << opt << "\t\t" << kStr << "\t" << desc << "\n";
}
void optPrint(std::ostringstream & os, std::string opt, bool keyval, std::vector<std::string> desc)
{
	std::string dstr = desc[0];
	for (int i = 1; i < desc.size(); ++i)
	{
		dstr += "\n\t\t\t" + desc[i];
	}
	optPrint(os,opt,keyval,dstr);
}
void helpMenu()
{
	std::ostringstream os;
	os << "=================================\n\n";
	os << "   Evenness Report Generator\n\n";
	os << "by Jack Fraser-Govil (2023)\n";
	os << "=================================\n\n";
	os << "Generate a typeset report, plots and csv files for the Poisson-evenness statistical measure, evaluated against the provided datasets\n\n";
	os << "Usage:\n\tlist-of-files | even [option] [option]....\n\n";
	os << "\tRun even on a piped list of files (such as from ls or cat). \n";
	os << "\tFiles can be relative locations. If no files provided, exits\n";
	os << "\teven must be run in a location where the following is true:\n";
	os << "\t\t-User has write access\n";
	os << "\t\t-pdflatex is available (unless -noreport used)\n";
	os << "\t\t-gnuplot is available (unless -noreport used and -plot unused)\n";
	os << "\tCoverage data must be space-delimited with:\n\t\t-col1 = the base coverage (integer) \n\t\t-col2 = the occurrence frequency (arbitrary normalisation)\n\n";
	os << "Options:\n";
	optPrint(os,"o",true,std::vector<std::string>{"Sets the name of the output files","Default value: report"});
	optPrint(os,"title",true,std::vector<std::string>{"Sets the title printed at the top of the report","Default value: `Evenness Report'"});
	optPrint(os,"a",true,(std::vector<std::string>){"Sets the location of the alias file","Gives human-readable names to the input files"});
	optPrint(os,"sort",false,std::vector<std::string>{"Turns on sorting on the output table by unevenness value","Default behaviour: use input file order"});
	optPrint(os,"noreport",false,"Turns off latex report generation");
	optPrint(os,"csv",false,"Turns on csv format output (in addition to other output)");
	optPrint(os,"plot",false,std::vector<std::string>{"If -noreport active, turns on figure plotting","Else, retains figures after report compilation"});
	optPrint(os,"v",false,"Toggles on verbose mode, including pdflatex output");
	optPrint(os,"open",false,"If report generation active, toggles on opening the output file using default viewer");
	std::cout << os.str() << std::endl;
}

void generateCrossRefs(std::string compileTitle,std::vector<std::string> fileList)
{
	std::string labelFile = compileTitle + "_annotations.txt";
	JSL::initialiseFile(labelFile);
	std::string refString = "#This is a generated alias file\n";
	refString += "#On each line, give the file an alias, the rerun the pipeline\n";
	refString += "#Valid alias cannot include special characters (including underscores), but can include space\n";
	refString += "#To omit a file from analysis, give it the alias REMOVE\n";
	JSL::writeStringToFile(labelFile,refString);

	JSL::writeVectorToFile(labelFile,fileList," \n",true);

	
}

int main(int argc, char**argv)
{
	JSL::Toggle help(false,"h",argc,argv);
	if (help)
	{
		helpMenu();
		exit(0);
	}

	JSL::Argument<std::string> CompileTitle("report","o",argc,argv);
	JSL::Argument<std::string> PrintTitle("__none__","title",argc,argv);
	JSL::Argument<std::string> crossRefFile("__none__","a",argc,argv);
	JSL::Toggle manualLabel(false,"label",argc,argv);
	JSL::Toggle sort(false,"sort",argc,argv);
	JSL::Toggle csv(false,"csv",argc,argv);
	JSL::Toggle plots(false,"plot",argc,argv);
	JSL::Toggle noreport(false,"noreport",argc,argv);
	JSL::Toggle loudCompile(false,"v",argc,argv);
	JSL::Toggle openReport(false,"open",argc,argv);
	bool keepFigs = plots;
	if (!noreport)
	{
		plots.Value = true;
	}

	std::vector<std::string> fileList;
	if (JSL::PipedInputFound())
	{
		
		forLineInPipedInput(
			fileList.push_back(PIPE_LINE);
		);
	}
	else
	{
		std::cout << "No piped input detected" << std::endl;
		exit(1);
	}

	std::vector<std::vector<std::string>> refs;
	if (manualLabel)
	{
		generateCrossRefs(CompileTitle,fileList);
		exit(0);
	}
	else
	{
		refs = getCrossRefs(crossRefFile);
	}
	

	std::vector<DataFile> fs;
	for (int i = 0; i < fileList.size(); ++i)
	{
		auto d = DataFile(fileList[i],refs);
		if (d.GoodFile)
		{
			fs.push_back(d);
		}
	}

	if (plots)
	{
		JSL::gnuplot gp;
		for (int i = 0; i < fs.size(); ++i)
		{
			fs[i].Plot(gp,i);
		}
		gp.SetFontSize(JSL::Fonts::Global,10);
		// gp.WindowSize(5,4);
		gp.SetXLog(true);
		gp.SetYLog(true);
		gp.SetTerminal("eps");
		gp.SetOutput((std::string)CompileTitle + "_plot.eps");
		gp.SetYLabel("Fractional Occurence");
		gp.SetXLabel("Base Coverage");
		gp.SetLegend(true);
		gp.SetFontSize(JSL::Fonts::Legend,8);
		gp.SetLegendLocation("outside");
		gp.SetLegendColumns(2);
		gp.Show();

		JSL::gnuplot gp2;
		gp2.SetFontSize(JSL::Fonts::Global,10);
		// gp2.WindowSize(5,4);
		for (int i = 0; i < fs.size(); ++i)
		{
			fs[i].TruncationPlot(gp2,i);
		}

		gp2.SetAxis(0);
		// gp2.SetXLog(true);
		// gp2.SetYLog(true);
		gp2.SetTerminal("eps");
		gp2.SetOutput((std::string)CompileTitle + "_plot_variance.eps");
		gp2.SetYLabel("Uneveness Score");
		gp2.SetXLabel("Truncation Factor (in units of Median Coverage)");
		gp2.SetLegend(true);
		gp2.SetFontSize(JSL::Fonts::Legend,8);
		gp2.SetLegendLocation("outside");
		gp.SetLegendColumns(2);
		gp2.Show();
	}
	
	if (!noreport)
	{
		PrintReport(CompileTitle,PrintTitle,sort,fs,loudCompile,openReport);
		std::string convCmd = "rm " + (std::string)CompileTitle + "*-eps-converted-to.pdf";
		system(convCmd.c_str());
		if (!keepFigs)
		{
			std::string cmd = "rm " + (std::string)CompileTitle + "*.eps";
			system(cmd.c_str());
		}
	}
	if (csv)
	{
		std::string dataFile = (std::string)CompileTitle + "_data.csv";
		JSL::initialiseFile(dataFile);
		JSL::writeStringToFile(dataFile,getTableString(fs));
	}
	return 0;
}