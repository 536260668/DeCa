//
// Created by 梦想家xixi on 2021/12/20.
//

#include "SAMBinaryTagAndValue.h"
#include "SAMUtils.h"
#include <stdexcept>

SAMBinaryTagAndValue::SAMBinaryTagAndValue(short tag, void *value, Void_Type voidType, int length) : tag(tag), value(value), type(voidType){
    if(value == nullptr) {
        throw std::invalid_argument("SAMBinaryTagAndValue value may not be null");
    } else if (!isAllowedAttributeValue(value, voidType)) {
        throw std::invalid_argument("Attribute type not supported.");
    } else {
        if((voidType != Uint8_t_Array_Type) && (voidType != Short_Array_Type) && (voidType != Int_Array_Type) && (voidType != Float_Array_Type))
            this->length = -1;
        else
            this->length = length;
    }
}

bool SAMBinaryTagAndValue::isAllowedAttributeValue(void* value, Void_Type voidType) {
    if((voidType != Uint8_Type) && (voidType != Short_Type) && (voidType != Integer_Type) && (voidType != String_Type) && (voidType != Character_Type)
    && (voidType != Float_Type) && (voidType != Uint8_t_Array_Type) && (voidType != Short_Array_Type) && (voidType != Int_Array_Type) && (voidType != Float_Array_Type)) {
        if((voidType != Long_Type)) {
            return false;
        } else {
            return SAMUtils::isValidUnsignedIntegerAttribute(*(long*)value) || *(long*)value >= -2147483648L && *(long*)value <= 2147483647L;
        }
    } else {
        return true;
    }
}

SAMBinaryTagAndValue *SAMBinaryTagAndValue::remove(SAMBinaryTagAndValue* root,short tag) {
    if(root->tag == tag){
        SAMBinaryTagAndValue* ret = root->next;
        delete root;
        return ret;
    }
    SAMBinaryTagAndValue* iter = root;
    while(iter->next != nullptr) {
        if(tag == iter->next->tag) {
            SAMBinaryTagAndValue* tmp = iter->next;
            iter->next = tmp->next;
            delete tmp;
            break;
        } else {
            iter = iter->next;
        }
    }
    return root;
}

SAMBinaryTagAndValue *SAMBinaryTagAndValue::insert(SAMBinaryTagAndValue *root, SAMBinaryTagAndValue *attr) {
    if(attr == nullptr)
        return root;
    else if (attr->next != nullptr) {
        throw std::invalid_argument("Can only insert single tag/value combinations.");
    } else if (attr->tag <= root->tag){
        attr->next = root;
        return attr;
    } else {
        SAMBinaryTagAndValue* iter = root;
        bool flag = false;
        while(iter->next != nullptr) {
            if(iter->next->tag > attr->tag) {
                SAMBinaryTagAndValue* tmp = iter->next;
                iter->next = attr;
                attr->next = tmp;
                flag = true;
                break;
            } else {
                iter = iter->next;
            }
        }
        if(!flag) {
            iter->next = attr;
            attr->next = nullptr;
        }
        return root;
    }
}

SAMBinaryTagAndValue *SAMBinaryTagAndValue::find(short tag) {
    SAMBinaryTagAndValue* iter = this;
    while(iter != nullptr && iter->tag <= tag) {
        if(iter->tag == tag) {
            return iter;
        } else {
            iter = iter->next;
        }
    }
    return nullptr;
}
