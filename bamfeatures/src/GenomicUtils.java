
public class GenomicUtils {

    public GenomicUtils() {

    }

    //reverseComplement
    public static String reverseComplement(String s){
        int i = 0;
        StringBuilder newString = new StringBuilder();
        //complement
        while(i < s.length()){
            if(s.charAt(i) == 'A'){
                newString.append('T');
            }
            if(s.charAt(i) == 'T'){
                newString.append('A');
            }
            if(s.charAt(i) == 'C'){
                newString.append('G');
            }
            if(s.charAt(i) == 'G'){
                newString.append('C');
            }
            i++;
        }
        //reverse
        newString.reverse();
        return newString.toString();
    }

    //mutate
    public static char mutate(char c){
        int rand = (int) (Math.random() * 2);

        if (c == 'A') {
            char[] n =  {'C', 'T', 'G'};
            return n[rand];
        }
        else if (c == 'C') {
            char[] n =  {'A', 'T', 'G'};
            return n[rand];
        }
        else if (c == 'G') {
            char[] n =  {'A', 'T', 'C'};
            return n[rand];
        }
        else if (c == 'T') {
            char[] n =  {'A', 'G', 'C'};
            return n[rand];
        }
        else{
            System.out.println("non valid char: " + c);
            return 'Q';
        }
    }


}
