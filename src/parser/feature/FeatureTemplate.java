package parser.feature;

import utils.Utils;

public class FeatureTemplate {
		
	public enum Arc {
	
		/**
		 * "H"	: head
		 * "M"	: modifier
		 * "B"	: in-between tokens
		 * 
		 * "P"	: pos tag
		 * "W"	: word form or lemma
		 * "EMB": word embedding (word vector)
		 * 
		 * "p": previous token
		 * "n": next token
		 *
		 */
		
		FEATURE_TEMPLATE_START,
	    
		/*************************************************
		 * Arc feature inspired by MST parser 
		 * ***********************************************/
		
	    // posL posIn posR
	    HP_BP_MP,			//CORE_POS_PC,	    
	    					//CORE_POS_XPC,
	    
	    // posL-1 posL posR posR+1
	    HPp_HP_MP_MPn,		//CORE_POS_PT0,
	    HP_MP_MPn,			//CORE_POS_PT1,
	    HPp_HP_MP,			//CORE_POS_PT2,
	    HPp_MP_MPn,			//CORE_POS_PT3,
	    HPp_HP_MPn,			//CORE_POS_PT4,
    
	    // posL posL+1 posR-1 posR
	    HP_HPn_MPp_MP,		//CORE_POS_APT0,
	    HP_MPp_MP,			//CORE_POS_APT1,
	    HP_HPn_MP,			//CORE_POS_APT2,
	    HPn_MPp_MP,			//CORE_POS_APT3,
	    HP_HPn_MPp,			//CORE_POS_APT4,
	    
	    // posL-1 posL posR-1 posR
	    // posL posL+1 posR posR+1
	    HPp_HP_MPp_MP,		//CORE_POS_BPT,
	    HP_HPn_MP_MPn,		//CORE_POS_CPT,

	    // unigram (form, lemma, pos, coarse_pos, morphology) 
	    CORE_HEAD_WORD,
	    CORE_HEAD_POS,
	    CORE_MOD_WORD,
	    CORE_MOD_POS,
	    CORE_HEAD_pWORD,
	    CORE_HEAD_nWORD,
	    CORE_MOD_pWORD,
	    CORE_MOD_nWORD,    
	    
	    // bigram  [word|lemma]-cross-[pos|cpos|mophlogy](-cross-distance)
	    HW_MW_HP_MP,			//CORE_BIGRAM_A,
	    MW_HP_MP,				//CORE_BIGRAM_B,
	    HW_HP_MP,				//CORE_BIGRAM_C,
	    MW_HP,					//CORE_BIGRAM_D,
	    HW_MP,					//CORE_BIGRAM_E,
	    HW_HP,					//CORE_BIGRAM_H,
	    MW_MP,					//CORE_BIGRAM_K,
	    HW_MW,					//CORE_BIGRAM_F,
	    HP_MP,					//CORE_BIGRAM_G,
	    
	    // label feature
	    CORE_LABEL_NTS1,
	    CORE_LABEL_NTH,
	    CORE_LABEL_NTI,
	    CORE_LABEL_NTIA,
	    CORE_LABEL_NTIB,
	    CORE_LABEL_NTIC,
	    CORE_LABEL_NTJ,

	    
		/*************************************************
		 * 2o feature  
		 * ***********************************************/

	    HP_SP_MP,
		HC_SC_MC,

		pHC_HC_SC_MC,
		HC_nHC_SC_MC,
		HC_pSC_SC_MC,
		HC_SC_nSC_MC,
		HC_SC_pMC_MC,
		HC_SC_MC_nMC,

		pHC_HL_SC_MC,
		HL_nHC_SC_MC,
		HL_pSC_SC_MC,
		HL_SC_nSC_MC,
		HL_SC_pMC_MC,
		HL_SC_MC_nMC,

		pHC_HC_SL_MC,
		HC_nHC_SL_MC,
		HC_pSC_SL_MC,
		HC_SL_nSC_MC,
		HC_SL_pMC_MC,
		HC_SL_MC_nMC,

		pHC_HC_SC_ML,
		HC_nHC_SC_ML,
		HC_pSC_SC_ML,
		HC_SC_nSC_ML,
		HC_SC_pMC_ML,
		HC_SC_ML_nMC,

		HC_MC_SC_pHC_pMC,
		HC_MC_SC_pHC_pSC,
		HC_MC_SC_pMC_pSC,
		HC_MC_SC_nHC_nMC,
		HC_MC_SC_nHC_nSC,
		HC_MC_SC_nMC_nSC,
		HC_MC_SC_pHC_nMC,
		HC_MC_SC_pHC_nSC,
		HC_MC_SC_pMC_nSC,
		HC_MC_SC_nHC_pMC,
		HC_MC_SC_nHC_pSC,
		HC_MC_SC_nMC_pSC,

		SP_MP,
		SW_MW,
		SW_MP,
		SP_MW,
		SC_MC,
		SL_ML,
		SL_MC,
		SC_ML,

		// head bigram
		H1P_H2P_M1P_M2P,
		H1P_H2P_M1P_M2P_DIR,
		H1C_H2C_M1C_M2C,
		H1C_H2C_M1C_M2C_DIR,

		// gp-p-c
		GP_HP_MP,
		GC_HC_MC,
		GL_HC_MC,
		GC_HL_MC,
		GC_HC_ML,

		GL_HL_MC,
		GL_HC_ML,
		GC_HL_ML,
		GL_HL_ML,

		GC_HC,
		GL_HC,
		GC_HL,
		GL_HL,

		GC_MC,	// this block only cross with dir flag
		GL_MC,
		GC_ML,
		GL_ML,
		HC_MC,
		HL_MC,
		HC_ML,
		HL_ML,

		pGC_GC_HC_MC,
		GC_nGC_HC_MC,
		GC_pHC_HC_MC,
		GC_HC_nHC_MC,
		GC_HC_pMC_MC,
		GC_HC_MC_nMC,

		pGC_GL_HC_MC,
		GL_nGC_HC_MC,
		GL_pHC_HC_MC,
		GL_HC_nHC_MC,
		GL_HC_pMC_MC,
		GL_HC_MC_nMC,

		pGC_GC_HL_MC,
		GC_nGC_HL_MC,
		GC_pHC_HL_MC,
		GC_HL_nHC_MC,
		GC_HL_pMC_MC,
		GC_HL_MC_nMC,

		pGC_GC_HC_ML,
		GC_nGC_HC_ML,
		GC_pHC_HC_ML,
		GC_HC_nHC_ML,
		GC_HC_pMC_ML,
		GC_HC_ML_nMC,

		GC_HC_MC_pGC_pHC,
		GC_HC_MC_pGC_pMC,
		GC_HC_MC_pHC_pMC,
		GC_HC_MC_nGC_nHC,
		GC_HC_MC_nGC_nMC,
		GC_HC_MC_nHC_nMC,
		GC_HC_MC_pGC_nHC,
		GC_HC_MC_pGC_nMC,
		GC_HC_MC_pHC_nMC,
		GC_HC_MC_nGC_pHC,
		GC_HC_MC_nGC_pMC,
		GC_HC_MC_nHC_pMC,

		// gp sibling
		GC_HC_MC_SC,
		GL_HC_MC_SC,
		GC_HL_MC_SC,
		GC_HC_ML_SC,
		GC_HC_MC_SL,

		// tri-sibling
		HC_PC_MC_NC,
		HL_PC_MC_NC,
		HC_PL_MC_NC,
		HC_PC_ML_NC,
		HC_PC_MC_NL,

		HC_PC_NC,
		PC_MC_NC,
		HL_PC_NC,
		HC_PL_NC,
		HC_PC_NL,
		PL_MC_NC,
		PC_ML_NC,
		PC_MC_NL,

		PC_NC,
		PL_NC,
		PC_NL,

		// ggpc
		GGC_GC_HC_MC,
		GGL_GC_HC_MC,
		GGC_GL_HC_MC,
		GGC_GC_HL_MC,
		GGC_GC_HC_ML,

		GGC_HC_MC,
		GGL_HC_MC,
		GGC_HL_MC,
		GGC_HC_ML,
		GGC_GC_MC,
		GGL_GC_MC,
		GGC_GL_MC,
		GGC_GC_ML,
		GGC_MC,
		GGL_MC,
		GGC_ML,
		GGL_ML,

		// psc
		HC_MC_CC_SC,
		HL_MC_CC_SC,
		HC_ML_CC_SC,
		HC_MC_CL_SC,
		HC_MC_CC_SL,

		HC_CC_SC,
		HL_CC_SC,
		HC_CL_SC,
		HC_CC_SL,

		// pp attachment
		PP_HC_MC,
		PP_HL_MC,
		PP_HC_ML,
		PP_HL_ML,

		PP_PL_HC_MC,
		PP_PL_HL_MC,
		PP_PL_HC_ML,
		PP_PL_HL_ML,

		// conjunction
		CC_CP_LP_RP,
		CC_CP_LC_RC,
		CC_CW_LP_RP,
		CC_CW_LC_RC,

		CC_LC_RC_FID,

		CC_CP_HC_AC,
		CC_CP_HL_AL,
		CC_CW_HC_AC,
		CC_CW_HL_AL,

		// PNX
		PNX_MW,
		PNX_HP_MW,

		// right branch
		RB,

		// child num
		CN_HP_NUM,
		CN_HL_NUM,
		CN_HP_LNUM_RNUM,
		CN_STR,

		// heavy
		HV_HP,
		HV_HC,

		// neighbor
		NB_HP_LC_RC,
		NB_HC_LC_RC,
		NB_HL_LC_RC,
		NB_GC_HC_LC_RC,
		NB_GC_HL_LC_RC,
		NB_GL_HC_LC_RC,

		// non-proj
		NP,
		NP_MC,
		NP_HC,
		NP_HL,
		NP_ML,
		NP_HC_MC,
		NP_HL_MC,
		NP_HC_ML,
		NP_HL_ML,

		/*************************************************
		 * word embedding feature  
		 * ***********************************************/
	    
	    HEAD_EMB,
	    MOD_EMB,
	    
	    
	    FEATURE_TEMPLATE_END;
		
		public static int numArcFeatBits = Utils.log2(FEATURE_TEMPLATE_END.ordinal());
	}

	public enum Word {

		FEATURE_TEMPLATE_START,
		
		/*************************************************
		 * Word features for matrix/tensor 
		 * ***********************************************/
		
		WORDFV_BIAS,
		
	    WORDFV_W0,
	    WORDFV_Wp,
	    WORDFV_Wn,
	    WORDFV_W0P0,
	    
	    WORDFV_F,
	    WORDFV_L0F,
	    WORDFV_P0F,
	    
	    WORDFV_VOICE,
	    WORDFV_L0VOICE,
	    
	    WORDFV_P0,
	    WORDFV_Pp,
	    WORDFV_Pn,
	    WORDFV_PpP0,
	    WORDFV_P0Pn,
	    WORDFV_PpP0Pn,
	    
	    WORDFV_EMB,
	    

	    
	    FEATURE_TEMPLATE_END;
	    
		public static int numWordFeatBits = Utils.log2(FEATURE_TEMPLATE_END.ordinal());
	}
	
	public enum Link {
		
		
		FEATURE_TEMPLATE_START,
		
		BIAS,
		
		PrW,
		PrP,
		ArW,
		ArP,
		
		PrW_ArW,
		PrP_ArP,
		PrW_ArP,
		PrP_ArW,

		
		PATH,
		PATH_PrW,
		PATH_PrP,
		PATH_ArW,
		PATH_ArP,
		PATH_PrP_ArP,
		PATH_PrW_ArP,
		PATH_ArW_PrP,
		PATH_PrW_ArW,
		
		PATHLEN,
		POSITION,
		
		// voice
		VOICE_DIR_PrP,
		VOICE_DIR_PrW,
		VOICE_DIR_PrP_ArP,
		
		// consecutive argument
		PP_A1P_A2P,
		PP_A1P_A2P_DIR,
		PW_A1P_A2P,
		PW_A1P_A2P_DIR,
		PP_A1W_A2P,
		PP_A1W_A2P_DIR,
		PP_A1P_A2W,
		PP_A1P_A2W_DIR,
		A1P_A2P,
		A1P_A2P_DIR,
		A1W_A2P,
		A1W_A2P_DIR,
		A1P_A2W,
		A1P_A2W_DIR,
		A1W_A2W,
		A1W_A2W_DIR,
		PW_A1W_A2P,
		PP_A1W_A2W,
		PW_A1P_A2W,
		
		// grandparent
		GPP_PP_AP,
		GPP_PP_AP_DIR,
		GPW_PP_AP,
		GPW_PP_AP_DIR,
		GPP_PW_AP,
		GPP_PW_AP_DIR,
		GPP_PP_AW,
		GPP_PP_AW_DIR,
		GPP_AP,
		GPP_AP_DIR,
		GPW_AP,
		GPW_AP_DIR,
		GPP_AW,
		GPP_AW_DIR,
		GPW_AW,
		GPW_AW_DIR,
		GPW_PW_AP,
		GPP_PW_AW,
		GPW_PP_AW,
		
		// co predicate
		PP_P2P_AP,
		PP_P2P_AP_DIR,
		PW_P2P_AP,
		PW_P2P_AP_DIR,
		PP_P2W_AP,
		PP_P2W_AP_DIR,
		PP_P2P_AW,
		PP_P2P_AW_DIR,
		PW_P2W_AP,
		PP_P2W_AW,
		PW_P2P_AW,
		
		// global
		VOICE_ARGSEQ,
		VOICE_LEMMA_ARGSEQ,
		VOICE_POSARGSEQ,
		VOICE_LEMMA_POSARGSEQ,
		VOICE_POSARGPREDSEQ,

		FEATURE_TEMPLATE_END;
		
		public static int numLinkFeatBits = Utils.log2(FEATURE_TEMPLATE_END.ordinal());
	}
	
	public enum Path {
		
		FEATURE_TEMPLATE_START,
		
		PATHFV_BIAS,
		
		PATHFV_PATHCODE,
		PATHFV_LABEL,
		PATHFV_DIS,
		PATHFV_DEP,
		PATHFV_POS,
		
		FEATURE_TEMPLATE_END;
		public static int numPathFeatBits = Utils.log2(FEATURE_TEMPLATE_END.ordinal());
	}
	
	public enum Context {
		
		FEATURE_TEMPLATE_START,
		
		CNTFV_BIAS,
		
		CNTFV_LABEL,
		
		FEATURE_TEMPLATE_END;
		public static int numContextFeatBits = Utils.log2(FEATURE_TEMPLATE_END.ordinal());
	}
}


