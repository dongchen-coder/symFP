//
//  idxAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/18/17.
//
//


#include "idxAnalysis.hpp"

using namespace llvm;

//#define IDX_DEBUG

namespace idxAnalysis {
    
    char IndexAnalysis::ID = 0;
    static RegisterPass<IndexAnalysis> X("idxAnalysis", "Array index analysis Pass", false, false);
    
    IndexAnalysis::IndexAnalysis() : FunctionPass(ID), arrayName(), arrayExpression() {}
    
    /* Pushing all instructions that related to the index expression of a load/store instruction into a stack (vector).
       This function will be called by function computeExpression() */
    std::vector<Instruction*> IndexAnalysis::instStackInit(Instruction* inst) {
        std::vector<Instruction*> instStack;
        instStack.push_back(inst);
        
        /* push all the instructions related to the target instruction into stack */
        unsigned long start = 0;
        while(instStack.size() != start) {
            
            if (instStack[start]->getNumOperands() == 1) {
                if (isa<Instruction>(instStack[start]->getOperand(0))) {
                    Instruction* pinst = dyn_cast<Instruction>(instStack[start]->getOperand(0));
                    if (std::find(instStack.begin(), instStack.end(), pinst) == instStack.end())     {
                        if (!isa<PHINode>(pinst) && !isa<AllocaInst>(pinst)) {
                            instStack.push_back(pinst);
                        }
                    }
                }
            }
            
            if (instStack[start]->getNumOperands() == 2) {
                if (isa<Instruction>(instStack[start]->getOperand(0))) {
                    Instruction* pinst = dyn_cast<Instruction>(instStack[start]->getOperand(0));
                    if (std::find(instStack.begin(), instStack.end(), pinst) == instStack.end()) {
                        if (!isa<PHINode>(pinst) && !isa<AllocaInst>(pinst)) {
                            instStack.push_back(pinst);
                        }
                    }
                }
                if (isa<Instruction>(instStack[start]->getOperand(1))) {
                    Instruction* pinst = dyn_cast<Instruction>(instStack[start]->getOperand(1));
                    if (std::find(instStack.begin(), instStack.end(), pinst) == instStack.end()) {
                        if (!isa<PHINode>(pinst) && !isa<AllocaInst>(pinst)) {
                            instStack.push_back(pinst);
                        }
                    }
                }
            }
            
            start++;
            if (start > 100) {
#ifdef IDX_DEBUG
                errs() << "ERROR: Expression should not be that long\n";
#endif
                break;
            }
        }
        return instStack;
    }
    
    
    /* Extracting array index expression from load/store instruction 
       1. Pushing all related instructions into stack by calling instStackInit()
       2. Scanning instruction stack to get a prefix expression 
       3. Convert prefix expression to infix expression */
    std::string IndexAnalysis::computeExpression(Instruction *inst, std::vector<std::string>& idvs) {
        
        std::vector<Instruction*> instStack = instStackInit(inst);
        
        if (instStack.size() == 0) {
            return "Warning: No expression";
        }
        
#ifdef IDX_DEBUG
        errs() << "============================\n";
        for (std::vector<Instruction*>::iterator it = instStack.begin(), eit = instStack.end(); it != eit; ++it) {
            (*it)->dump();
        }
        errs() << "============================\n";
#endif
        
        /* reverse iterating stack to get expression */
        std::map<Instruction*, std::string> instExprs;
        std::vector<std::string> opExprs;
        for (std::vector<Instruction*>::reverse_iterator it = instStack.rbegin(), eit = instStack.rend(); it != eit; ++it) {
            Instruction* instStackEntry = *it;
            
            if (isa<AllocaInst>(instStackEntry) || isa<LoadInst>(instStackEntry)) {
                continue;
            }
            
            /* Extract expressions for operands */
            opExprs.clear();
            for (unsigned int i = 0; i < instStackEntry->getNumOperands(); i++) {
                if (isa<Argument>(instStackEntry->getOperand(i))) {
                    opExprs.push_back(instStackEntry->getOperand(i)->getName().str());
                } else if (isa<Constant>(instStackEntry->getOperand(i))) {
                    Constant* cons = dyn_cast<Constant>(instStackEntry->getOperand(i));
                    opExprs.push_back(cons->getUniqueInteger().toString(10, true));
                } else if (isa<LoadInst>(instStackEntry->getOperand(i))) {
                    LoadInst* ld = dyn_cast<LoadInst>(instStackEntry->getOperand(i));
                    opExprs.push_back(ld->getOperand(0)->getName().str());
                    idvs.push_back(ld->getOperand(0)->getName().str());
                } else if (isa<PHINode>(instStackEntry->getOperand(i))) {
                    if (instStackEntry->getOperand(i)->hasName()) {
                        opExprs.push_back(instStackEntry->getOperand(i)->getName().str());
                    } else {
#ifdef IDX_DEBUG
                        errs() << "Error: no name for PHInode\n";
#endif
                    }
                } else if (isa<Instruction>(instStackEntry->getOperand(i))) {
                    Instruction* instTmp = dyn_cast<Instruction>(instStackEntry->getOperand(i));
                    if (instExprs.find(instTmp) != instExprs.end()) {
                        opExprs.push_back(instExprs[instTmp]);
                    } else {
                        opExprs.push_back("");
#ifdef IDX_DEBUG
                        errs() << "Error: unsolved previous instructions\n";
#endif
                    }
                } else {
#ifdef IDX_DEBUG
                    instStackEntry->getType()->dump();
                    errs() << "Warning: unknown operand type for instruction in InstStack\n";
#endif
                }
            }
            
            /* Assamble current instruction */
            std::string exprForCurrentInst = "";
            if (instStackEntry->isCast()) {
                exprForCurrentInst = opExprs[0];
            } else if (instStackEntry->isBinaryOp()) {
#ifdef IDX_DEBUG
                errs() << instStackEntry->getNumOperands() << " " << opExprs.size() << "\n";
#endif
                BinaryOperator* BinaryOpTmp = dyn_cast<BinaryOperator>(instStackEntry);
                exprForCurrentInst += "(";
                exprForCurrentInst += opExprs[0];
                switch (BinaryOpTmp->getOpcode()) {
                    case llvm::Instruction::Add:
                        exprForCurrentInst += " + ";
                        break;
                    case llvm::Instruction::FAdd:
                        exprForCurrentInst += " + ";
                        break;
                    case llvm::Instruction::Sub:
                        exprForCurrentInst += " - ";
                        break;
                    case llvm::Instruction::FSub:
                        exprForCurrentInst += " - ";
                        break;
                    case llvm::Instruction::Mul:
                        exprForCurrentInst += " * ";
                        break;
                    case llvm::Instruction::FMul:
                        exprForCurrentInst += " * ";
                        break;
                    case llvm::Instruction::FDiv:
                        exprForCurrentInst += " / ";
                        break;
                    case llvm::Instruction::UDiv:
                        exprForCurrentInst += " / ";
                        break;
                    case llvm::Instruction::SDiv:
                        exprForCurrentInst += " / ";
                        break;
                    default:
                        break;
                }
                exprForCurrentInst += opExprs[1];
                exprForCurrentInst += ")";
            } else if (instStackEntry->isShift()) {
                //TODO
            } else if (isa<CallInst>(instStackEntry)) {
                //TODO
            } else {
#ifdef IDX_DEBUG
                errs() << "Error: unregconized instruction in InstStack\n";
#endif
                return "";
            }
#ifdef IDX_DEBUG
            instStackEntry->dump();
            errs() << exprForCurrentInst << "\n";
#endif
            
            /* Push the assambled instruction to map */
            instExprs[instStackEntry] = exprForCurrentInst;
        }
        
        return instExprs[inst];
    }
    
    /* Extracting array name from load/store instruction */
    std::string IndexAnalysis::getArrayName(GetElementPtrInst* inst) {
        std::string nameTmp = "";
        
        if (isa<LoadInst>(inst->getPointerOperand())) {
            // this array is a read, cost = 1
            /* Load instruction */
            LoadInst* ldTmp = dyn_cast<LoadInst>(inst->getPointerOperand());
            nameTmp = ldTmp->getOperand(0)->getName().str();
        } else {
            // this array is a write, cost = 2
            /* Store instruction */
            nameTmp = inst->getPointerOperand()->getName().str();
        }
        
        return nameTmp;
    }
    
    /* Scan all IR instructions in function to find load/stores accessing arrays */
    void IndexAnalysis::findAllArrayAccesses(Function &F) {
        /* Using instruction iterator to interating through all the IR instructions in the functions */
        for (inst_iterator it = inst_begin(F), eit = inst_end(F); it != eit; ++it) {
            int cost = 0;
            if (isa<LoadInst>(*it)) {
                std::vector<std::string> idvTmp;
                if (isa<GetElementPtrInst>(it->getOperand(0))) {
                    GetElementPtrInst* gepTmp = dyn_cast<GetElementPtrInst>(it->getOperand(0));
                    // this is a load array reference
                    cost += 1;

                    /* Problem: how to figure out which operand is index for GEP? */
                    int OperandToAnalysis = 1;
                    
                    if (gepTmp->getNumOperands() < 2) {
#ifdef IDX_DEBUG
                        errs() << "ERROR: GEP does not have enough operands\n";
#endif
                    } else {
                        Instruction* accessInst = dyn_cast<Instruction>(&(*it));
                        std::string nameTmp = getArrayName(gepTmp);
                        
                        if (isa<Instruction>(gepTmp->getOperand(OperandToAnalysis))) {
                            Instruction* indexInst = dyn_cast<Instruction>(gepTmp->getOperand(OperandToAnalysis));
                            std::string idxExpr = computeExpression(indexInst, idvTmp);
                            // this is a indirect access expression
                            if (!idxExpr.empty()/*expr is not empty*/) 
                                cost += 1;
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = idxExpr;
                            arrayAccessVariable[accessInst] = idvTmp;
                        } else if (isa<Argument>(gepTmp->getOperand(OperandToAnalysis))) {
                            // this is a indirect access expression
                            cost += 1;
                            Argument* argTmp = dyn_cast<Argument>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = argTmp->getName().str();
                        } else if (isa<Constant>(gepTmp->getOperand(OperandToAnalysis))) { 
                            // this is a constant access expression
                            cost += 0;
                            Constant* constTmp = dyn_cast<Constant>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = constTmp->getUniqueInteger().toString(10, true);
                        } else {
#ifdef IDX_DEBUG
                            errs() << "ERROR: other load address type\n";
#endif
                        }
                    }
                    BCReference[gepTmp] = cost;
                }
            }
            
            if (isa<StoreInst>(*it)) {
                std::vector<std::string> idvTmp;
                if (isa<GetElementPtrInst>(it->getOperand(1))) {
                    GetElementPtrInst* gepTmp = dyn_cast<GetElementPtrInst>(it->getOperand(1));
                    // this is a load array reference
                    cost += 2;
                    
                    /* Problem: how to figure out which operand is index for GEP? */
                    int OperandToAnalysis = 1;
                    
                    if (gepTmp->getNumOperands() < 2) {
#ifdef IDX_DEBUG
                        errs() << "ERROR: GEP does not have enough operands\n";
#endif
                    } else {
                        Instruction* accessInst = dyn_cast<Instruction>(&(*it));
                        std::string nameTmp = getArrayName(gepTmp);
                        if (isa<Instruction>(gepTmp->getOperand(OperandToAnalysis))) {
                            Instruction* indexInst = dyn_cast<Instruction>(gepTmp->getOperand(OperandToAnalysis));
                            std::string idxExpr = computeExpression(indexInst, idvTmp);
                            // this is a indirect access expression
                            if (!idxExpr.empty()/*expr is not empty*/) 
                                cost += 1;
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = idxExpr;
                            arrayAccessVariable[accessInst] = idvTmp;
                        } else if (isa<Argument>(gepTmp->getOperand(OperandToAnalysis))) {
                            // this is a indirect access expression
                            cost += 1;
                            Argument* argTmp = dyn_cast<Argument>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = argTmp->getName().str();
                        } else if (isa<Constant>(gepTmp->getOperand(OperandToAnalysis))) {
                             // this is a constant access expression
                            cost += 0;
                            Constant* constTmp = dyn_cast<Constant>(gepTmp->getOperand(OperandToAnalysis));
                            arrayName[accessInst] = nameTmp;
                            arrayExpression[accessInst] = constTmp->getUniqueInteger().toString(10, true);
                        } else {
#ifdef IDX_DEBUG
                            errs() << "ERROR: other store address type\n";
#endif
                        }
                    }
                    BCReference[gepTmp] = cost;
                }
            }
        }
        return;
    }
    
    /* Dump (array name, index expression) pairs for all references inside the function.
       Two hash tables are used to store array names and index expressions, the keys are load/store instruction pointers */
    void IndexAnalysis::dumpAllInfo() {
        errs() << "Array index info: Total number of references: " << arrayName.size()  << "\n";
        for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
            errs() << it->second << " " << arrayExpression[it->first] << "\n";
        }
        double average_cost = 0.0;
        errs() << "BC Array cost info: Total number of arrays: " << BCArrayList.size() << "\n";
        for(std::map<std::string, int>::iterator it = BCArrayList.begin(); it != BCArrayList.end(); ++it) {
            average_cost += it->second;
            errs() << it->first << " " << it->second << "\n";
        }
        errs() << "BC Average cost: " << average_cost / BCArrayList.size() << "\n";
        return;
    }
    
    /* Analysis main function: information are dumped as comments at the beginning of generated static sampling code */
    bool IndexAnalysis::runOnFunction(Function &F) {
        errs() << "\n /* Start to analysis array index\n";
        
        findAllArrayAccesses(F);

        // merge cost value from reference based to array based
        for(std::map<Instruction*, int>::iterator it = BCReference.begin(); it != BCReference.end(); ++it) {
            if (isa<GetElementPtrInst>(it->first)) {
                GetElementPtrInst* gepTmp = dyn_cast<GetElementPtrInst>(it->first);
                std::string name = getArrayName(gepTmp);
                if (BCArrayList.find(name) != BCArrayList.end()) 
                    BCArrayList[name] += it->second;
                else
                    BCArrayList[name] = it->second;
            }
        }

        BCArrayRatio = BCArrayRatio / BCArrayList.size();
        
        dumpAllInfo();
        
        errs() << "\n Finish to analysis array index */ \n";
        
        return false;
    }
    
}


