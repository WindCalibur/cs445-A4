import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Scanner;
import java.util.Arrays;

// Global alignment algorithm using affine gap penalty - general algorithm

public class A4Q3 {
	String mSeqA;
	String mSeqB;
	
	static String S1;
	static String S2;
	
	// ...
	double[][] M,Gx,Gy;
	int[][] Mp,Gxp,Gyp;
	double[] q = new double[256];
	static double[][] emission;
	int mScore;
	String mAlignmentSeqA = "";
	String mAlignmentSeqB = "";
	double tow,eps,delta;
	final static Charset ENCODING = StandardCharsets.US_ASCII;
	final static String delims = "[ ]+";
	final static double NEGINF = Double.NEGATIVE_INFINITY / 2.0;

	static void init_sim(String matrix_name) throws IOException {
		emission = new double[256][256];
		for (int i=0; i<256; i++)
			for (int j=0; j<256; j++)
				emission[i][j] = 0;
		Path path = Paths.get(matrix_name);
		try (Scanner scanner = new Scanner(path,ENCODING.name())) {
			String line;
			do {
				line = scanner.nextLine();
			} while (line.charAt(0)=='#');
			String[] letters = line.split(delims);
			int map[];
			map = new int[21];
			for (int j=1; j<=20; j++) {
				map[j] = letters[j].charAt(0);
				emission[0][map[j]] = 1;
			}
			for (int i=0; i<20; i++) {
				String[] scores = scanner.nextLine().split(delims);
				int offset = 0;
				while (scores[offset].isEmpty())
					offset++;
				int row = scores[offset].charAt(0);
				for (int j=1; j<=20; j++) {
					emission[row][map[j]] = Double.parseDouble(scores[offset+j]);
				}
			}
			// copy C to U
			emission[0]['U'] = 1;
			for (int j=1; j<=20; j++) {
				emission['U'][map[j]] = emission['C'][map[j]];
				emission[map[j]]['U'] = emission[map[j]]['C'];
			}
		}
	}

	static String readName(Scanner scanner) {
		if (scanner.hasNextLine()) {
			String[] parts = scanner.nextLine().split(delims);
			String[] subparts = parts[0].split("\\|");
			return subparts[2];
		}
		return null;
	}

	static String readSequence(Scanner scanner) {
		String output = new String();
		while (scanner.hasNextLine() && scanner.findInLine(">") == null)
			output += scanner.nextLine();
		return output;
	}

	static boolean checkSequence(String seq) {
		for (int i=0; i<seq.length(); i++)
			if (emission[0][seq.charAt(i)] == 0) {
				System.out.println("Illegal character: "+seq.charAt(i));
				return false;
			}
		return true;
	}

	void init(String seqA, String seqB) {
		mSeqA = seqA;
		mSeqB = seqB;
		mAlignmentSeqA = "";
		mAlignmentSeqB = "";
		M = new double[mSeqA.length() + 1][mSeqB.length() + 1];
		Gx = new double[mSeqA.length() + 1][mSeqB.length() + 1];
		Gy = new double[mSeqA.length() + 1][mSeqB.length() + 1];
		Mp = new int[mSeqA.length() + 1][mSeqB.length() + 1];
		Gxp = new int[mSeqA.length() + 1][mSeqB.length() + 1];
		Gyp = new int[mSeqA.length() + 1][mSeqB.length() + 1];
		M[0][0] = 0;
		Gx[0][0] = NEGINF;
		Gy[0][0] = NEGINF;
	}

	double process() {
//		System.out.println(mSeqA);
//		System.out.println(mSeqB);
		for (int i = 0; i <= mSeqA.length(); i++) {
			for (int j = 0; j <= mSeqB.length(); j++) {

				if ((i == 0) && (j == 0)) {
					continue;
				}
				
				// Prev, i-1,j-1
				double m = NEGINF;
				double x = NEGINF;
				double y = NEGINF;
				M[i][j] = NEGINF;
				if ((i != 0) && (j != 0)) {
					m = M[i-1][j-1] + Math.log(1-2*delta-tow);
					x = Gx[i-1][j-1]+ Math.log(1-eps-tow);
					y = Gy[i-1][j-1] + Math.log(1-eps-tow);
					M[i][j] = Math.log(emission[mSeqA.charAt(i-1)][mSeqB.charAt(j-1)]);
				}
				double temp = m;
				Mp[i][j] = 0;
				if (x > m) {
					temp = x;
					Mp[i][j] = 1;
				}
				if (y > temp) {
					temp = y;
					Mp[i][j] = 2;
				}
				M[i][j] = M[i][j] + temp;
				
				// Prev, i-1,j
				m = NEGINF;
				x = NEGINF;
				Gx[i][j] = NEGINF;
				if (i != 0) {
					m = M[i-1][j] + Math.log(delta);
					x = Gx[i-1][j] + Math.log(eps);
					Gx[i][j] = Math.log(q[mSeqA.charAt(i-1)]);
				}				
				Gxp[i][j] = 0;
				temp = m;
				if (x > m) {
					temp = x;
					Gxp[i][j] = 1;
				}
				
				Gx[i][j] = Gx[i][j] + temp;
				
				// Prev, i-1,j
				m = NEGINF;
				y = NEGINF;
				Gy[i][j] = NEGINF;
				if (j != 0) {
					m = M[i][j-1] + Math.log(delta);
					y = Gy[i][j-1]+ Math.log(eps);
					Gy[i][j] = Math.log(q[mSeqB.charAt(j-1)]);
				}
				
				Gyp[i][j] = 0;
				temp = m;
				if (y > m) {
					temp = y;
					Gyp[i][j] = 2;
				}
				
				Gy[i][j] = Gy[i][j] + temp;
				
//				System.out.println("M is " + " (" + i + "," + j + "): " +  M[i][j]);
//				System.out.println("Gx is " + " (" + i + "," + j + "): " +   Gx[i][j]);
//				System.out.println("Gy is " + " (" + i + "," + j + "): " +   Gy[i][j]);
			}
		}

		double score = M[mSeqA.length()][mSeqB.length()];
		if (Gx[mSeqA.length()][mSeqB.length()] > score)
			score = Gx[mSeqA.length()][mSeqB.length()];
		if (Gy[mSeqA.length()][mSeqB.length()] > score)
			score = Gy[mSeqA.length()][mSeqB.length()];
		
		return score + Math.log(tow);
	}

	
	void traceback(String seqA, String seqB) {
		// changes S1, S2
		
		int currentI = seqA.length();
		int currentJ = seqB.length();
		S1 = "";
		S2 = "";

		int state = 0;
		double score = M[seqA.length()][seqB.length()];
		if (Gx[seqA.length()][seqB.length()] > score) {
			score = Gx[seqA.length()][seqB.length()];
			state = 1;
		}

		if (Gy[seqA.length()][seqB.length()] > score) {
			score = Gy[seqA.length()][seqB.length()];
			state = 2;
		}
		// 0 = M, 1 = Gx, 2 = Gy
		while (true) {
			if ((currentI == 0) && (currentJ == 0)) {
				break;
			}
			int pointer = Mp[currentI][currentJ];
			if (state == 1) {
				pointer = Gxp[currentI][currentJ];
			} else if (state == 2) {
				pointer = Gyp[currentI][currentJ];
			}

			if (state == 0) {
				S1 = seqA.charAt(currentI-1) + S1;
				S2 = seqB.charAt(currentJ-1) + S2;
				currentI--;
				currentJ--;
			} else if (state == 2) {
				S1 = '-' + S1;
				S2 = seqB.charAt(currentJ-1) + S2;
				currentJ--;
			} else if (state == 1) {
				S1 = seqA.charAt(currentI-1) + S1;
				S2 = '-' + S2;
				currentI--;
			}
			state = pointer;
			
		}
		return;
	}
	 
	void printMatrix() {
		System.out.println("M =");
		for (int i = 0; i < mSeqA.length() + 1; i++) {
			for (int j = 0; j < mSeqB.length() + 1; j++) {
				System.out.print(String.format("%4d ", M[i][j]));
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Gx =");
		for (int i = 0; i < mSeqA.length() + 1; i++) {
			for (int j = 0; j < mSeqB.length() + 1; j++) {
				System.out.print(String.format("%4d ", Gx[i][j]));
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Gy =");
		for (int i = 0; i < mSeqA.length() + 1; i++) {
			for (int j = 0; j < mSeqB.length() + 1; j++) {
				System.out.print(String.format("%4d ", Gy[i][j]));
			}
			System.out.println();
		}
		System.out.println();
	}

	void printScoreAndAlignments() {
		System.out.println("Score: " + mScore);
		System.out.println("Sequence A: " + mAlignmentSeqB);
		System.out.println("Sequence B: " + mAlignmentSeqA);
		System.out.println();
	}

	A4Q3(double delta, double eps, double tow) {
		this.delta = delta;
		this.eps = eps;
		this.tow = tow;
		q['A'] = 5.99e-02;
		q['R'] = 5.60e-02;
		q['N'] = 4.82e-02;
		q['D'] = 5.22e-02;
		q['C'] = 3.31e-02;
		q['U'] = 3.31e-02;
		q['Q'] = 5.02e-02;
		q['E'] = 5.15e-02;
		q['G'] = 5.97e-02;
		q['H'] = 3.62e-02;
		q['I'] = 5.96e-02;
		q['L'] = 6.08e-02;
		q['K'] = 5.33e-02;
		q['M'] = 5.15e-02;
		q['F'] = 5.02e-02;
		q['P'] = 5.07e-02;
		q['S'] = 5.47e-02;
		q['T'] = 5.46e-02;
		q['W'] = 1.93e-02;
		q['Y'] = 4.01e-02;
		q['V'] = 5.83e-02;
		
	}
	
	public static void main(String [] args) {
		
		class Record implements Comparable<Record> {
			public int index;
			public String name;
			public String sequence;
			public double score;
			int length;
			double[][] M,Gx,Gy;
			String firstSequence;
			String secondSequence;

			@Override
			public int compareTo(Record r1) {
		        if (this.score < r1.score) return 1;
		        if (this.score > r1.score) return -1;
		        return 0;
		    }
		};
		
		String query;
		Record[] record = new Record[1017];
		try {
			init_sim("emission.txt");
//			Path path = Paths.get("unknown.txt");
//			Scanner scanner = new Scanner(path,ENCODING.name());
//			query = readSequence(scanner);
//			if (!checkSequence(query))
//				System.out.println("Query sequence contains illegal character");
			Path path = Paths.get("2017-01-16 uniprot.fasta");
			Scanner scanner = new Scanner(path,ENCODING.name());
			int maxlength = 0, maxlengthi = -1;
			
			for (int i=0; i<record.length; i++) {
				record[i] = new Record();
				record[i].index = i;
				record[i].name = readName(scanner);
				record[i].sequence = readSequence(scanner);
				if (!checkSequence(record[i].sequence))
					System.out.println(record[i].sequence+" contains illegal character");
				else
					if (record[i].sequence.length() > maxlength) {
						maxlength = record[i].sequence.length();
						maxlengthi = i;
					}
			}
			query = record[1000].sequence;
			System.out.println(maxlengthi+" L="+maxlength);
		} catch (IOException e) {
			System.out.println("Cannot find file "+e.toString());
			return;
		}

		A4Q3 nw = new A4Q3(0.08,0.35,0.002);
		

		for (int i=0; i<50; i++) {
			nw.init(record[i].sequence, record[59].sequence);
			record[i].score = nw.process();
			nw.traceback(record[i].sequence,record[59].sequence);
			record[i].firstSequence = S1;
			record[i].secondSequence = S2;
		}
	
		Arrays.sort(record);
		
		System.out.println("Sample Input and Output \n");
		int offset = (record.length-50);
		for (int i=offset; i<3+offset; i++) {
			System.out.println("Index=" + offset+i + " Name=" + record[i].name + " ln Pr=" + record[i].score);
			System.out.println(record[i].firstSequence);
			System.out.println(record[i].secondSequence + "\n");
		}
		
		// Real data
		
		for (int i=0; i<1000; i++) {
			nw.init(record[i].sequence, query);
			record[i].score = nw.process();
			nw.traceback(record[i].sequence,query);
			record[i].firstSequence = S1;
			record[i].secondSequence = S2;
		}
	
		Arrays.sort(record);
		offset = (record.length-1000);
		System.out.println("Real Input and Output \n");
		for (int i=offset; i<offset+3; i++) {
			System.out.println("Index=" + offset+i + " Name=" + record[i].name + " ln Pr=" + record[i].score);
			System.out.println(record[i].firstSequence);
			System.out.println(record[i].secondSequence + "\n");
		}
		
	}
}
