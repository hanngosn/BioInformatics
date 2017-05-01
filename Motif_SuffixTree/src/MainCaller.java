/*
Created on April 28, 2017

@author: Hannie
*/

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;




public class MainCaller {
	
	public static void main(String[] args) {
		Tree in = new Tree();
		
		String[] words = new String[] {"ACGTCTAGTA","TAGTAAAA", "TTAAGGTCTAGTA","TCTAGTAGA"};
        for (int i = 0; i < words.length; ++i) 
            in.put(words[i], i);
        
               
//        find most common motif 
        int k = 4; 	//kmer length
        int n = 3;	//number of motifs
        		
        ArrayList<String> ls_motif= new ArrayList<String>();

		
        for (int i = 0; i < words.length; ++i) 
            for (String s : getSubstrings(words[i])) {
            	if (s.length() == k && !Arrays.asList(ls_motif).contains(s)){
            		System.out.println(in.search(s).size());	
            		ls_motif.add(s);
            	}
            }
        System.out.println(ls_motif);
        
	}    
        
        
	public static Set<String> getSubstrings(String str) {
        Set<String> ret = new HashSet<String>();
        
        // compute all substrings
        for (int len = 1; len <= str.length(); ++len) {
            for (int start = 0; start + len <= str.length(); ++start) {
                String itstr = str.substring(start, start + len);
                ret.add(itstr);
            }
        }
        System.out.println(ret);
        return ret;
    }
	
}
