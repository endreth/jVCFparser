package model;

import java.util.Map;

public class VCF {

    private Map vcfMap;
    private int popID;

    public VCF() {
        this.vcfMap = vcfMap;
        this.popID = 1;
    }

    public Map getVcfMap() {
        return vcfMap;
    }

    public void setVcfMap(Map vcfMap) {
        this.vcfMap = vcfMap;
    }
}
